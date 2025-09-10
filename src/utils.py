import os
from pathlib import Path

import numpy as np
import pandas as pd
from MSCI.Grouping_MS1.Grouping_mw_irt import (
    make_data_compatible,
    find_combinations_kdtree,
    within_ppm,
    within_tolerance,
)
from MSCI.Preprocessing.Koina import PeptideProcessor
from MSCI.Similarity.spectral_angle_similarity import joinPeaks, nspectraangle
from matchms.importing import load_from_msp
from scipy.spatial import cKDTree


def read_pool_file(fname):
    df = pd.read_csv(fname, usecols=["Pool_name", "Sequence"])
    df = df[~df["Sequence"].str.contains("-")]
    df = df[df["Sequence"].str.len() <= 30]
    df = df.drop_duplicates("Sequence").reset_index(drop=True)
    return df


def read_pool_files(fnames):
    return (
        pd.concat([read_pool_file(f) for f in fnames])
        .drop_duplicates("Sequence")
        .reset_index(drop=True)
    )


def read_pool_dir(dir_path, ext=".txt"):
    all_files = [
        os.path.join(dir_path, f) for f in os.listdir(dir_path) if f.endswith(ext)
    ]
    return read_pool_files(all_files)


def write_db(db, out_path, col_name="Sequence"):
    os.makedirs(os.path.dirname(out_path), exist_ok=True)
    db[col_name].to_csv(out_path, index=False, header=False)


def predict_spectra(
    input_file: str,
    collision_energy: int = 30,
    charge: int = 2,
    model_intensity: str = "Prosit_2020_intensity_HCD",
    model_irt: str = "Prosit_2019_irt",
):
    pred_file = f"{Path(input_file).with_suffix('')}.msp"

    if os.path.exists(pred_file):
        spectra = list(load_from_msp(pred_file))
        return spectra, pred_file

    processor = PeptideProcessor(
        input_file=input_file,
        collision_energy=collision_energy,
        charge=charge,
        model_intensity=model_intensity,
        model_irt=model_irt,
    )

    processor.process(pred_file)

    if Path(pred_file).stat().st_size == 0:
        raise RuntimeError("Generating spectrum predictions failed")

    spectra = list(load_from_msp(pred_file))
    return spectra, pred_file


def process_spectra_pairs(
    chunk, spectra, mz_irt_df, tolerance=0, ppm=0, m=0, n=0.5, progress_queue=None
):
    results = []

    for index_pair in chunk:
        i, j = index_pair

        x = spectra[i]
        y = spectra[j]

        x_df = pd.DataFrame({"mz": x.peaks.mz, "intensities": x.peaks.intensities})
        y_df = pd.DataFrame({"mz": y.peaks.mz, "intensities": y.peaks.intensities})

        matcher = joinPeaks(tolerance=tolerance, ppm=ppm)
        x_matched, y_matched = matcher.match(x_df, y_df)

        angle = nspectraangle(x_matched, y_matched, m=m, n=n)

        # Extract the relevant information for the given index pair
        results.append(
            {
                "index1": i,
                "index2": j,
                "peptide 1": mz_irt_df.loc[i, "Name"],
                "peptide 2": mz_irt_df.loc[j, "Name"],
                "m/z  1": mz_irt_df.loc[i, "MW"],
                "m/z 2": mz_irt_df.loc[j, "MW"],
                "iRT 1": mz_irt_df.loc[i, "iRT"],
                "iRT 2": mz_irt_df.loc[j, "iRT"],
                "similarity_score": angle,
            }
        )

        if progress_queue is not None:
            progress_queue.put(1)  # report progress

    return pd.DataFrame(results)


def process_search_db_pairs(
    chunk,
    spectra_search,
    spectra_db,
    mz_irt_search_df,
    mz_irt_db_df,
    tolerance=0,
    ppm=0,
    m=0,
    n=0.5,
    progress_queue=None,
):
    results = []

    for search_i, db_i in chunk:
        x = spectra_search[search_i]
        y = spectra_db[db_i]

        x_df = pd.DataFrame({"mz": x.peaks.mz, "intensities": x.peaks.intensities})
        y_df = pd.DataFrame({"mz": y.peaks.mz, "intensities": y.peaks.intensities})

        matcher = joinPeaks(tolerance=tolerance, ppm=ppm)
        x_matched, y_matched = matcher.match(x_df, y_df)

        angle = nspectraangle(x_matched, y_matched, m=m, n=n)

        # Extract the relevant information for the given index pair
        results.append(
            {
                "index_search": search_i,
                "index_db": db_i,
                "peptide_search": mz_irt_search_df.loc[search_i, "Name"],
                "peptide_db": mz_irt_db_df.loc[db_i, "Name"],
                "m/z_search": mz_irt_search_df.loc[search_i, "MW"],
                "m/z_db": mz_irt_db_df.loc[db_i, "MW"],
                "iRT_search": mz_irt_search_df.loc[search_i, "iRT"],
                "iRT_db": mz_irt_db_df.loc[db_i, "iRT"],
                "similarity_score": angle,
            }
        )

        if progress_queue is not None:
            progress_queue.put(1)  # report progress

    return pd.DataFrame(results)


def process_peptide_combinations(mz_irt_df, tolerance1, tolerance2, use_ppm=True):
    compatible_data = make_data_compatible(mz_irt_df)
    result_tolerance = find_combinations_kdtree(
        compatible_data, tolerance1, tolerance2, use_ppm
    )
    unique_result_tolerance = list({tuple(sorted(pair)) for pair in result_tolerance})

    # Create a DataFrame for the results
    results = []
    for (index1, mw1, irt1), (index2, mw2, irt2) in unique_result_tolerance:
        results.append(
            {
                "index1": index1,
                "index2": index2,
                "peptide 1": mz_irt_df.loc[index1, "Name"],
                "peptide 2": mz_irt_df.loc[index2, "Name"],
                "m/z  1": mz_irt_df.loc[index1, "MW"],
                "m/z 2": mz_irt_df.loc[index2, "MW"],
                "iRT 1": mz_irt_df.loc[index1, "iRT"],
                "iRT 2": mz_irt_df.loc[index2, "iRT"],
            }
        )

    results_df = pd.DataFrame(results)
    return results_df


def find_combinations_kdtree_cross(
    search_data, query_data, tolerance1, tolerance2, use_ppm=True
):
    """
    For each point in search_data, find similar points in query_data based on tolerance values.

    Parameters:
        search_data: list of tuples (index, mw, irt)
        query_data: list of tuples (index, mw, irt)
        tolerance1, tolerance2: tolerances for mw and irt
        use_ppm: if True, tolerance1 is interpreted as ppm

    Returns:
        List of tuples: [(search_point, query_point), ...]
    """
    valid_combinations = []

    # Convert datasets to numpy arrays for k-d tree
    search_array = np.array([(mw, irt) for index, mw, irt in search_data])
    query_array = np.array([(mw, irt) for index, mw, irt in query_data])

    # Build k-d tree for the query dataset
    tree = cKDTree(query_array)

    if use_ppm:
        ppm_tolerance1_sq = (tolerance1 / 1e6) ** 2
        tolerance2_sq = tolerance2**2

    # Iterate through each point in the search dataset
    for i, point in enumerate(search_array):
        if use_ppm:
            radius = np.sqrt(ppm_tolerance1_sq * point[0] ** 2 + tolerance2_sq)
        else:
            radius = np.sqrt(tolerance1**2 + tolerance2**2)

        # Find all points in the query dataset within the radius
        indices = tree.query_ball_point(point, radius)

        for j in indices:
            pair = (search_data[i], query_data[j])
            if use_ppm:
                if within_ppm(pair, tolerance1, tolerance2):
                    valid_combinations.append(pair)
            else:
                if within_tolerance(pair, tolerance1, tolerance2):
                    valid_combinations.append(pair)

    return valid_combinations


def process_peptide_db_combinations(
    mz_irt_search_df, mz_irt_db_df, tolerance1, tolerance2, use_ppm=True
):
    compatible_search_data = make_data_compatible(mz_irt_search_df)
    compatible_db_data = make_data_compatible(mz_irt_db_df)

    result_tolerance = find_combinations_kdtree_cross(
        compatible_search_data, compatible_db_data, tolerance1, tolerance2, use_ppm
    )
    unique_result_tolerance = list({tuple(sorted(pair)) for pair in result_tolerance})

    # Create a DataFrame for the results
    results = []
    for (index1, mw1, irt1), (index2, mw2, irt2) in unique_result_tolerance:
        results.append(
            {
                "index_search": index1,
                "index_db": index2,
                "peptide_search": mz_irt_search_df.loc[index1, "Name"],
                "peptide_db": mz_irt_db_df.loc[index2, "Name"],
                "m/z_search": mz_irt_search_df.loc[index1, "MW"],
                "m/z_db": mz_irt_db_df.loc[index2, "MW"],
                "iRT_search": mz_irt_search_df.loc[index1, "iRT"],
                "iRT_db": mz_irt_db_df.loc[index2, "iRT"],
            }
        )

    results_df = pd.DataFrame(results)
    return results_df
