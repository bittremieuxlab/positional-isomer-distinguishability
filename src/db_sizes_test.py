import os
import time
from collections import defaultdict

import pandas as pd
from matplotlib import pyplot as plt

from generate_random_dbs import generate_random_dbs
from msci_runner import find_indistinguishable_peptides
from generate_positional_isomers import random_positional_isomer
from utils import write_db, read_pool_dir


def time_test():
    sizes = [1e3, 2e3, 4e3, 6e3, 8e3, 1e4, 2e4]
    times = []
    for s in sizes:
        t = time.time()
        print(f"Testing a db size: {s}")
        db = pd.read_csv("data/generated_positional_isomers.csv", header=None)
        db = db.sample(n=int(s))
        pep_file = f"random_dbs/pos_iso_{s}.txt"
        write_db(db, out_path=pep_file, col_name=0)
        find_indistinguishable_peptides(pep_file)
        print(f"Size {s} took {time.time() - t:.2f} seconds")
        times.append(time.time() - t)

    print(sizes, times)
    plt.plot(
        sizes, times, marker="o", linestyle="-"
    )  # optional: add markers for clarity
    plt.xlabel("DB size")
    plt.ylabel("Runtime (s)")
    plt.grid(True)
    plt.show()


def replace_with_positional_isomers(df, frac=0.5):
    to_replace = df.sample(frac=frac)
    to_keep = df[~df.index.isin(to_replace.index)]

    # Take the number of positional isomers to generate from to_keep
    to_add = to_keep.sample(n=len(to_replace), replace=True)
    to_add["Sequence"] = to_add["Sequence"].apply(random_positional_isomer)
    return pd.concat([to_keep, to_add])


def db_size_diversity(dir_path="data", sizes=None, diversities=None, out_dir=""):
    random_dbs = generate_random_dbs(dir_path=dir_path, sizes=sizes)
    os.makedirs(out_dir, exist_ok=True)
    results = defaultdict(dict)
    for db in random_dbs:
        for diversity in diversities:
            out_file = os.path.join(out_dir, f"{len(db)}_{diversity}.txt")
            if not os.path.isfile(out_file):
                div_db = replace_with_positional_isomers(db, frac=diversity)
                write_db(div_db, out_file)

            indistinguisable_peptides_csv = os.path.join(
                out_dir, f"{len(db)}_{diversity}.csv"
            )
            if not os.path.isfile(indistinguisable_peptides_csv):
                find_indistinguishable_peptides(out_file)

            results[len(db)][diversity] = pd.read_csv(indistinguisable_peptides_csv)
    return results


def plot_size_diversity_dict(data):
    """
    Plot values from nested dict structure.

    Parameters
    ----------
    data : dict
        {size: {diversity: DataFrame}}
    """

    records = []
    for size, diversities in data.items():
        for diversity, df in diversities.items():
            val = len(df[df["similarity_score"] >= 0.7]) * 2 / size
            records.append({"size": size, "diversity": diversity, "value": val})

    plot_df = pd.DataFrame(records)

    plt.figure(figsize=(8, 5))
    for diversity, group in plot_df.groupby("diversity"):
        plt.plot(
            group["size"], group["value"], marker="o", label=f"Diversity={diversity}"
        )

    plt.xlabel("Size")
    plt.ylabel("Value")
    plt.title("Value vs Size for Different Diversities")
    plt.legend()
    plt.grid(True)
    plt.show()


if __name__ == "__main__":
    # time_test()
    results = db_size_diversity(
        sizes=[5e4, 1e5, 1.5e5],
        diversities=[0.0, 0.25, 0.5, 0.75],
        out_dir="db_size_diversity",
    )
    plot_size_diversity_dict(results)
