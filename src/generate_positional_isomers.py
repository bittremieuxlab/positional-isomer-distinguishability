import pandas as pd
import random

from argparse import ArgumentParser


def parse_args():
    parser = ArgumentParser(description="Find indistinguishable peptides using MSCI")
    parser.add_argument("--input", "-i", required=True, help="Input file with peptides")
    parser.add_argument(
        "--output", "-o", required=True, help="Output file for generated isomers"
    )
    return parser.parse_args()


def random_positional_isomer(peptide: str) -> str:
    """Generate a random positional isomer by swapping two random amino acides."""
    if len(peptide) < 2:
        return peptide  # No swap possible

    # pick two distinct positions
    i, j = random.sample(range(len(peptide)), 2)

    # swap characters
    s_list = list(peptide)
    s_list[i], s_list[j] = s_list[j], s_list[i]

    return "".join(s_list)


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
    args = parse_args()
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
