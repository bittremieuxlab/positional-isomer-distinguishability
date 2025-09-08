from pathlib import Path
import pandas as pd

from MSCI.Grouping_MS1.Grouping_mw_irt import process_peptide_combinations
from MSCI.Preprocessing.Koina import PeptideProcessor
from MSCI.Preprocessing.Parsing import read_msp_file
from MSCI.Similarity.spectral_angle_similarity import process_spectra_pairs
from matchms.importing import load_from_msp

from argparse import ArgumentParser

parser = ArgumentParser(description="Find indistinguishable peptides using MSCI")
parser.add_argument("--input", "-i", required=True, help="Input file with peptides")
parser.add_argument(
    "--output", "-o", required=True, help="Output file for generated isomers"
)
args = parser.parse_args()


def generate_positional_isomers(peptide: str) -> list[str]:
    """
    Generate all positional isomers of a given peptide by moving each amino acid
    to every possible position in the peptide.
    :param peptide: The original peptide sequence.
    :return: A list of positional isomers.
    """
    isomers = set()
    length = len(peptide)
    for i in range(length):
        for j in range(length):
            if i != j:
                # Move amino acid from position i to position j
                new_pep = list(peptide)
                aa = new_pep.pop(i)
                new_pep.insert(j, aa)
                isomers.add("".join(new_pep))
    return list(isomers)


if __name__ == "__main__":
    input_file = args.input
    output_file = args.output
    df = pd.read_csv(input_file, header=None, names=["Peptide"])
    isomers = []
    for peptide in df["Peptide"]:
        isomers.extend(generate_positional_isomers(peptide))
    df_isomers = pd.DataFrame(isomers, columns=["Peptide"])
    df_isomers = df_isomers.drop_duplicates().reset_index(drop=True)
    df_isomers.to_csv(
        output_file,
        index=False,
        header=False,
    )
