import pandas as pd
from itertools import combinations
from collections import defaultdict
from concurrent.futures import ThreadPoolExecutor, as_completed
from tqdm import tqdm
from src.utils import read_pool_dir

import random
import os


def are_positional_isomers(seq1: str, seq2: str) -> bool:
    """
    Check if two peptides are positional isomers:
    - Must be same length
    - Must differ in exactly two positions
    - The amino acids at those positions must be swapped
    """
    diffs = [(a, b) for a, b in zip(seq1, seq2) if a != b]
    # Must differ in exactly two positions, and be swapped
    return len(diffs) == 2 and diffs[0] == diffs[1][::-1]


def group_by_composition(peptides: list[str]) -> dict[str, list[str]]:
    """
    Group peptides by their amino acid composition.
    Key = sorted string of amino acids (multiset signature).
    """
    groups = defaultdict(list)
    for pep in peptides:
        signature = "".join(sorted(pep))
        groups[signature].append(pep)
    return groups


def find_isomers_in_group(peptides: list[str], pbar=None) -> list[tuple[str, str]]:
    isomer_pairs = []
    for pep1, pep2 in combinations(peptides, 2):
        if are_positional_isomers(pep1, pep2):
            isomer_pairs.append((pep1, pep2))
        if pbar:
            pbar.update(1)
    return isomer_pairs


def find_positional_isomers(
    peptides: list[str], max_workers: int = 4
) -> list[tuple[str, str]]:
    """
    Find all pairs of peptides in the list that are positional isomers.
    Uses grouping by composition to reduce comparisons.
    """
    isomer_pairs = []
    groups = group_by_composition(peptides)
    # Count total comparisons
    total_comparisons = sum(len(list(combinations(g, 2))) for g in groups.values())
    with tqdm(total=total_comparisons, desc="Checking pairs") as pbar:
        with ThreadPoolExecutor(max_workers=max_workers) as executor:
            futures = [
                executor.submit(find_isomers_in_group, g, pbar)
                for g in groups.values()
                if len(g) > 1
            ]
            for future in as_completed(futures):
                isomer_pairs.extend(future.result())
    return isomer_pairs


def main(data_folder: str = "data"):
    df = read_pool_dir(data_folder)
    isomers = find_positional_isomers(df["Sequence"].unique().tolist(), max_workers=8)
    df_isomers = pd.DataFrame(isomers, columns=["Peptide1", "Peptide2"])
    # Put both columns under one column and drop duplicates
    df_all_isomers = (
        pd.concat([df_isomers["Peptide1"], df_isomers["Peptide2"]])
        .drop_duplicates()
        .reset_index(drop=True)
    )
    df_all_isomers.to_csv(
        "positional_isomers_msci_input.csv",
        index=False,
        header=False,
    )


if __name__ == "__main__":
    main()
