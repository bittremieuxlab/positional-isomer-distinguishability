import logging
from functools import partial
from multiprocessing import cpu_count, Pool
from pathlib import Path

import pandas as pd
from MSCI.Preprocessing.Koina import PeptideProcessor
from MSCI.Preprocessing.Parsing import read_msp_file
from matchms.importing import load_from_msp
from utils import process_spectra_pairs, process_peptide_combinations

from argparse import ArgumentParser

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


def parallelize(func, data, n_chunks=None, **kwargs):
    if n_chunks is None:
        n_chunks = cpu_count()

    # Split into chunks
    chunk_size = (len(data) + n_chunks - 1) // n_chunks
    chunks = [data[i : i + chunk_size] for i in range(0, len(data), chunk_size)]

    # Use partial to pass extra args
    func_with_args = partial(func, **kwargs)

    with Pool(n_chunks) as pool:
        dfs = pool.map(func_with_args, chunks)

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
    processor = PeptideProcessor(
        input_file=input_file,
        collision_energy=collision_energy,
        charge=charge,
        model_intensity=model_intensity,
        model_irt=model_irt,
    )

    pred_file = f"{Path(input_file).with_suffix('')}.msp"

    if output_file is None:
        if input_file.endswith(".csv"):
            output_file = f"{Path(input_file).with_suffix('')}.csv.csv"
        else:
            output_file = f"{Path(input_file).with_suffix('')}.csv"

    processor.process(pred_file)

    if Path(pred_file).stat().st_size == 0:
        raise RuntimeError("Generating spectrum predictions failed")

    print("Got spectrum predictions")
    spectra = list(load_from_msp(pred_file))
    mz_irt_df = read_msp_file(pred_file)
    groups_df = process_peptide_combinations(
        mz_irt_df, mz_tolerance, irt_tolerance, use_ppm=False
    )
    print("Got peptide groups")
    groups_df.columns = groups_df.columns.str.replace(" ", "")
    index_array = groups_df[["index1", "index2"]].values.astype(int)
    result = parallelize(
        process_spectra_pairs,
        index_array,
        n_chunks=n_chunks,
        spectra=spectra,
        mz_irt_df=mz_irt_df,
        tolerance=peak_tolerance,
        ppm=peak_ppm,
    )
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
