import logging
from argparse import ArgumentParser
from functools import partial
from multiprocessing import cpu_count, Pool, Manager
from pathlib import Path

import pandas as pd
from MSCI.Preprocessing.Parsing import read_msp_file
from tqdm import tqdm

from utils import (
    predict_spectra,
    process_search_db_pairs,
    process_spectra_pairs,
    process_peptide_combinations,
)

logging.getLogger("matchms").setLevel(logging.ERROR)


def parse_args():
    parser = ArgumentParser(description="Find indistinguishable peptides using MSCI")
    parser.add_argument("--input", "-i", required=True, help="Input file with peptides")
    parser.add_argument("--output", "-o", help="Output file")
    parser.add_argument(
        "--collision_energy", "-ce", type=int, default=30, help="Collision energy"
    )
    parser.add_argument("--charge", "-c", type=int, default=2, help="Charge")
    parser.add_argument(
        "--model_intensity",
        "-mi",
        type=str,
        default="Prosit_2020_intensity_HCD",
        help="Model for intensity prediction",
    )
    parser.add_argument(
        "--model_irt",
        "-mirt",
        type=str,
        default="Prosit_2019_irt",
        help="Model for iRT prediction",
    )
    parser.add_argument(
        "--mz_tolerance",
        "-mzt",
        type=float,
        default=1,
        help="m/z tolerance for grouping peptides",
    )
    parser.add_argument(
        "--irt_tolerance",
        "-irtt",
        type=float,
        default=5,
        help="iRT tolerance for grouping peptides",
    )
    parser.add_argument(
        "--peak_tolerance",
        "-pt",
        type=float,
        default=0,
        help="Tolerance for peak matching in similarity calculation",
    )
    parser.add_argument(
        "--peak_ppm",
        "-ppm",
        type=float,
        default=10,
        help="PPM tolerance for peak matching in similarity calculation",
    )
    parser.add_argument(
        "--n-chunks",
        "-nc",
        type=int,
        default=None,
        help="Number of chunks to process in parallel, default (None) uses your cpu count.",
    )
    return parser.parse_args()


def parallel_process_spectra_pairs(
    spectra_pairs, spectra, mz_irt_df, tolerance=0, ppm=0, m=0, n=0.5, n_chunks=None
):
    """
    spectra_pairs: list of (i,j) index tuples
    spectra: list of spectra objects
    mz_irt_df: DataFrame with peptide info
    """

    if n_chunks is None:
        n_chunks = cpu_count()

    # Split the pairs into roughly equal chunks
    chunk_size = (len(spectra_pairs) + n_chunks - 1) // n_chunks
    chunks = [
        spectra_pairs[i : i + chunk_size]
        for i in range(0, len(spectra_pairs), chunk_size)
    ]

    with Manager() as manager:
        progress_queue = manager.Queue()
        func = partial(
            process_spectra_pairs,
            spectra=spectra,
            mz_irt_df=mz_irt_df,
            tolerance=tolerance,
            ppm=ppm,
            m=m,
            n=n,
            progress_queue=progress_queue,
        )

        with Pool(n_chunks) as pool:
            # Launch jobs
            results_async = pool.map_async(func, chunks)

            # Show progress bar
            with tqdm(
                total=len(spectra_pairs), desc="Processing spectra pairs"
            ) as pbar:
                processed = 0
                while processed < len(spectra_pairs):
                    progress_queue.get()
                    processed += 1
                    pbar.update(1)

            # Collect results
            dfs = results_async.get()

    return pd.concat(dfs, ignore_index=True)


def parallel_process_search_db_pairs(
    spectra_pairs,
    spectra_search,
    spectra_db,
    mz_irt_search_df,
    mz_irt_db_df,
    tolerance=0,
    ppm=0,
    m=0,
    n=0.5,
    n_chunks=None,
):
    if n_chunks is None:
        n_chunks = cpu_count()

    # Split the pairs into roughly equal chunks
    chunk_size = (len(spectra_pairs) + n_chunks - 1) // n_chunks
    chunks = [
        spectra_pairs[i : i + chunk_size]
        for i in range(0, len(spectra_pairs), chunk_size)
    ]

    with Manager() as manager:
        progress_queue = manager.Queue()
        func = partial(
            process_search_db_pairs,
            spectra_search=spectra_search,
            spectra_db=spectra_db,
            mz_irt_search_df=mz_irt_search_df,
            mz_irt_db_df=mz_irt_db_df,
            tolerance=tolerance,
            ppm=ppm,
            m=m,
            n=n,
            progress_queue=progress_queue,
        )

        with Pool(n_chunks) as pool:
            # Launch jobs
            results_async = pool.map_async(func, chunks)

            # Show progress bar
            with tqdm(
                total=len(spectra_pairs), desc="Processing spectra pairs"
            ) as pbar:
                processed = 0
                while processed < len(spectra_pairs):
                    progress_queue.get()
                    processed += 1
                    pbar.update(1)

            # Collect results
            dfs = results_async.get()

    return pd.concat(dfs, ignore_index=True)


def find_indistinguishable_peptides(
    input_file: str,
    collision_energy: int = 30,
    charge: int = 2,
    model_intensity: str = "Prosit_2020_intensity_HCD",
    model_irt: str = "Prosit_2019_irt",
    mz_tolerance: int | float = 1,
    irt_tolerance: int | float = 5,
    peak_tolerance: int | float = 0,
    peak_ppm: int | float = 10,
    output_file: str = None,
    n_chunks: int = None,
):
    """
    From a given file with peptides find those that could be indistinguishable

    :param input_file: text file with a single peptide on each line
    :param collision_energy: used for intensity prediction
    :param charge: used for intensity prediction
    :param model_intensity: Koina model name used for intensity prediction
    :param model_irt: Koina model name used for intensity prediction
    :param mz_tolerance: tolerance for mz difference for peptides to be considered in similarity calculation
    :param irt_tolerance: tolerance for irt difference for peptides to be considered in similarity calculation
    :param peak_tolerance: tolerance used when matching peaks
    :param peak_ppm: also used for tolerance when matching peaks
    :param output_file: path to output file or, if None, reuse the path and filename of the input file with different extension
    :param n_chunks: number of chunks to process in parallel, default (None) uses your cpu count.

    :return:
    """
    spectra, pred_file = predict_spectra(
        input_file=input_file,
        collision_energy=collision_energy,
        charge=charge,
        model_intensity=model_intensity,
        model_irt=model_irt,
    )

    mz_irt_df = read_msp_file(pred_file)
    groups_df = process_peptide_combinations(
        mz_irt_df, mz_tolerance, irt_tolerance, use_ppm=False
    )

    groups_df.columns = groups_df.columns.str.replace(" ", "")
    index_array = groups_df[["index1", "index2"]].values.astype(int)
    result = parallel_process_spectra_pairs(
        index_array,
        n_chunks=n_chunks,
        spectra=spectra,
        mz_irt_df=mz_irt_df,
        tolerance=peak_tolerance,
        ppm=peak_ppm,
    )

    if output_file is None:
        if input_file.endswith(".csv"):
            output_file = f"{Path(input_file).with_suffix('')}.csv.csv"
        else:
            output_file = f"{Path(input_file).with_suffix('')}.csv"

    result.to_csv(output_file, index=False)


if __name__ == "__main__":
    args = parse_args()
    find_indistinguishable_peptides(
        input_file=args.input,
        collision_energy=args.collision_energy,
        charge=args.charge,
        model_intensity=args.model_intensity,
        model_irt=args.model_irt,
        mz_tolerance=args.mz_tolerance,
        irt_tolerance=args.irt_tolerance,
        peak_tolerance=args.peak_tolerance,
        peak_ppm=args.peak_ppm,
        output_file=args.output,
        n_chunks=args.n_chunks,
    )
