from pathlib import Path

from MSCI.Grouping_MS1.Grouping_mw_irt import process_peptide_combinations
from MSCI.Preprocessing.Koina import PeptideProcessor
from MSCI.Preprocessing.Parsing import read_msp_file
from MSCI.Similarity.spectral_angle_similarity import process_spectra_pairs
from matchms.importing import load_from_msp


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
    min_similarity: float = 0.7,
    output_file: str = None,
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
    :param min_similarity: Filter on final similarity score
    :param output_file: path to output file or, if None, reuse the path and filename of the input file with different extension
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

    spectra = list(load_from_msp(pred_file))
    mz_irt_df = read_msp_file(pred_file)
    groups_df = process_peptide_combinations(
        mz_irt_df, mz_tolerance, irt_tolerance, use_ppm=False
    )
    groups_df.columns = groups_df.columns.str.replace(" ", "")
    index_array = groups_df[["index1", "index2"]].values.astype(int)
    result = process_spectra_pairs(
        index_array, spectra, mz_irt_df, tolerance=peak_tolerance, ppm=peak_ppm
    )
    result = result[result["similarity_score"] >= min_similarity]
    result.to_csv(output_file, index=False)


if __name__ == "__main__":
    find_indistinguishable_peptides("random_dbs/1000.txt")
