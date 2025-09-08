import os

import pandas as pd
from MSCI.Grouping_MS1.Grouping_mw_irt import (
    make_data_compatible,
    find_combinations_kdtree,
)
from MSCI.Similarity.spectral_angle_similarity import joinPeaks, nspectraangle


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


def process_spectra_pairs(chunk, spectra, mz_irt_df, tolerance=0, ppm=0, m=0, n=0.5):
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
