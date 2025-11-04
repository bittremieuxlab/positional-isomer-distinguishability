import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
import scipy.stats as stats
import ast
from matplotlib.ticker import FuncFormatter


def is_positional_isomer(peptide1: str, peptide2: str) -> bool:
    return peptide1 != peptide2 and sorted(peptide1) == sorted(peptide2)


def process(val):
    if val is None:
        return None
    aa_str = "".join(ast.literal_eval(str(val)))
    return "".join(sorted(aa_str))


def label_positional_isomers(df, pep1="peptide 1", pep2="peptide 2") -> pd.DataFrame:
    is_positional_isomer_list = []
    positional_isomer_distance_list = []
    positional_aa_list = []
    first_difference_position_list = []
    second_difference_position_list = []
    length_list = []
    for peptide1, peptide2 in zip(df[pep1], df[pep2]):
        diffs = [(a, b) for a, b in zip(peptide1, peptide2) if a != b]
        is_positional_isomer_list.append(len(diffs) == 2 and diffs[0] == diffs[1][::-1])
        if diffs == []:
            positional_aa_list.append(None)
        else:
            positional_aa_list.append(diffs[0])
        # Must differ in exactly two positions, and be swapped
        if (len(diffs) == 2 and diffs[0] == diffs[1][::-1]) == True:
            length = len(peptide1)
            diff_dist = [i for i in range(length) if peptide1[i] != peptide2[i]]
            positional_isomer_distance_list.append(diff_dist[1] - diff_dist[0])
            first_difference_position_list.append(diff_dist[0])
            second_difference_position_list.append(diff_dist[1])
            length_list.append(length)
        else:
            positional_isomer_distance_list.append(None)
            first_difference_position_list.append(None)
            second_difference_position_list.append(None)
            length_list.append(None)
    df["is_positional_isomer"] = is_positional_isomer_list
    df["positional_amino_acids"] = positional_aa_list
    df["is_positional_isomer_distance"] = positional_isomer_distance_list
    df["first_difference_position"] = first_difference_position_list
    df["second_difference_position"] = second_difference_position_list
    df["peptide_length"] = length_list
    df["pep1_IL"] = df[pep1].str.replace("I", "L")
    df["pep2_IL"] = df[pep2].str.replace("I", "L")
    df["aa_content_1_IL"] = df["pep1_IL"].apply(lambda x: "".join(sorted(x)))
    df["aa_content_2_IL"] = df["pep2_IL"].apply(lambda x: "".join(sorted(x)))
    df["positional_aa"] = df["positional_amino_acids"].apply(process)
    return df


def plot_isomer_distance(df, plot_path: str = None):
    plt.figure(figsize=(9, 6))
    sns.boxplot(
        data=df,
        x="is_positional_isomer_distance",
        y="similarity_score",
    )
    plt.title("Similarity Score by Positional Isomer Distance")
    plt.xlabel("Positional Isomer Distance")
    plt.ylabel("Similarity Score")
    plt.ylim(0, 1)
    plt.grid(axis="y")
    sns.despine()
    plt.gca().xaxis.set_major_formatter(FuncFormatter(lambda x, _: f"{int(x)}"))
    plt.savefig(plot_path + "/isomer_distance.png", dpi=300)


def plot_isomer_position(df):
    plt.figure(figsize=(10, 6))
    sns.boxplot(
        data=df,
        x="first_difference_position",
        y="similarity_score",
    )
    plt.title("Similarity Score by Positional Isomer Position")
    plt.xlabel("Positional Isomer Position")
    plt.ylabel("Similarity Score")
    plt.ylim(0, 1)
    plt.grid(axis="y")
    plt.show()


def plot_isomer_aa(df):
    # Order df based on median similarity score per positional_aa
    medians = (
        df.groupby("positional_aa")["similarity_score"]
        .median()
        .sort_values(ascending=False)
    )
    plt.figure(figsize=(10, 6))
    sns.boxplot(
        data=df,
        order=medians.index,
        x="positional_aa",
        y="similarity_score",
    )
    plt.title("Similarity Score by Positional Isomer Amino Acids")
    plt.xlabel("Positional Isomer Amino Acids")
    plt.ylabel("Similarity Score")
    plt.ylim(0, 1)
    plt.xticks(rotation=90)
    plt.grid(axis="y")
    plt.show()


def report_isomer_stats(df):
    print(f"Total pairs: {len(df)}")
    print(f"Pairs with similarity == 1: {len(df[df['similarity_score'] == 1])}")
    print(
        f"Pairs with the exact same peptide: {len(df[df['peptide_search'] == df['peptide_db']])}"
    )
    print(
        f"Pairs with distinct peptides and a similarity > 0.7: {len(df[~(df['peptide_search'] == df['peptide_db']) & (df['similarity_score'] > 0.7)])}"
    )
    print(f"Positional isomers: {df['is_positional_isomer'].sum()}")
    print(
        f"Positional isomers with similarity > 0.7: {len(df[(df['is_positional_isomer'] == True) & (df['similarity_score'] > 0.7)])}"
    )
    print(
        f"I/L substitutions with similarity > 0.7: {len(df[~(df['peptide_search'] == df['peptide_db']) & (df['pep1_IL'] == df['pep2_IL']) & (df['similarity_score'] > 0.7)])}"
    )
    print(
        f"Same AA composition but not positional isomers with similarity > 0.7: {len(df[(df['pep1_IL'] != df['pep2_IL']) & (df['is_positional_isomer'] == False) & (df['similarity_score'] > 0.7) & (df['aa_content_1_IL']==df['aa_content_2_IL'])])}"
    )
    print(
        f"Positional isomer counts: {df[(df['similarity_score'] > 0.7)& (df['is_positional_isomer'] == True)]['positional_aa'].value_counts()}"
    )
    peptide_match_list = df[df["peptide_search"] == df["peptide_db"]][
        "peptide_search"
    ].drop_duplicates()
    found_peptides = df[(df["similarity_score"] > 0.7)][
        "peptide_search"
    ].drop_duplicates()
    distinct_found_peptides = found_peptides[~found_peptides.isin(peptide_match_list)]
    print(
        f"Number of peptides found that did not have an exact match: {len(distinct_found_peptides)}"
    )


def compare_two_outputs(
    path1, path2, rt_cutoff=5, similarity_cutoff=0.7, pep1="peptide 1", pep2="peptide 2"
):
    df1 = pd.read_csv(path1)
    df2 = pd.read_csv(path2)
    df1 = label_positional_isomers(df1)
    df2 = label_positional_isomers(df2)
    df1["RT_diff"] = np.abs(df1["iRT 1"] - df1["iRT 2"])
    df2["RT_diff"] = np.abs(df2["iRT 1"] - df2["iRT 2"])
    df1 = df1[df1["RT_diff"] < rt_cutoff]
    df2 = df2[df2["RT_diff"] < rt_cutoff]
    print(f"Path 1: {len(df1)} pairs after RT cutoff")
    print(f"Path 2: {len(df2)} pairs after RT cutoff")
    df1_high_sim = df1[df1["similarity_score"] > similarity_cutoff]
    df2_high_sim = df2[df2["similarity_score"] > similarity_cutoff]
    print(f"Path 1: {len(df1_high_sim)} pairs after similarity cutoff")
    print(f"Path 2: {len(df2_high_sim)} pairs after similarity cutoff")
    # Merge on peptide 1 and peptide 2
    merged = df1_high_sim.merge(
        df2_high_sim,
        on=[pep1, pep2],
        suffixes=("_path1", "_path2"),
    )
    print(f"Merged: {len(merged)} pairs in both")
    # Plot similarity scores against each other
    plt.figure(figsize=(8, 8))
    sns.scatterplot(
        data=merged,
        x="similarity_score_path1",
        y="similarity_score_path2",
        hue="is_positional_isomer_path1",
        alpha=0.7,
    )
    plt.plot([0, 1], [0, 1], "k--", label="y=x")
    plt.xlim(0, 1)
    plt.ylim(0, 1)
    plt.xlabel("Similarity Score Path 1")
    plt.ylabel("Similarity Score Path 2")
    plt.title("Comparison of Similarity Scores")
    plt.legend()
    plt.grid()
    plt.show()
    # List of peptides only in one of the two
    only_in_path1 = df1_high_sim[
        ~df1_high_sim.set_index([pep1, pep2]).index.isin(
            df2_high_sim.set_index([pep1, pep2]).index
        )
    ]
    only_in_path2 = df2_high_sim[
        ~df2_high_sim.set_index([pep1, pep2]).index.isin(
            df1_high_sim.set_index([pep1, pep2]).index
        )
    ]
    print(f"Only in Path 1: {len(only_in_path1)} pairs")
    print(f"Only in Path 2: {len(only_in_path2)} pairs")


def plot_aa_grid(df):
    # Split into separate AAs
    df[["aa1", "aa2"]] = df["positional_aa"].apply(lambda x: pd.Series(list(x)))
    # Median similarity per unique pair
    median_sim = df.groupby(["aa1", "aa2"])["similarity_score"].median().reset_index()
    # Pivot into symmetric matrix
    grid = median_sim.pivot(index="aa1", columns="aa2", values="similarity_score")
    # Ensure full symmetric matrix (fill both triangles)
    grid = grid.combine_first(grid.T)
    plt.figure(figsize=(8, 6))
    cmap = plt.cm.viridis
    im = plt.imshow(grid, cmap=cmap, interpolation="nearest")
    plt.colorbar(im, label="Median Similarity")
    # Annotate only visible values
    # for (i, j), val in np.ndenumerate(grid.values):
    #     plt.text(j, i, f"{val:.2f}", ha="center", va="center", color="w")
    plt.xticks(range(len(grid.columns)), grid.columns)
    plt.yticks(range(len(grid.index)), grid.index)
    plt.title("AA Swap Median Similarity")
    plt.show()


def plot_aa_grid_bottom_triangle(df, plot_path):
    aa_groups = {
        "Nonpolar": ["A", "V", "L", "I", "M", "F", "W", "P", "G"],
        "Polar (uncharged)": ["S", "T", "C", "Y", "N", "Q"],
        "Positive (basic)": ["K", "R", "H"],
        "Negative (acidic)": ["D", "E"],
    }
    ordered_aas = [aa for group in aa_groups.values() for aa in group]
    # Split into separate AAs
    df[["aa1", "aa2"]] = df["positional_aa"].apply(lambda x: pd.Series(list(x)))
    # Count per unique pair
    counts = df.groupby(["aa1", "aa2"]).size().reset_index(name="count")
    grid = counts.pivot(index="aa1", columns="aa2", values="count")
    grid = grid.combine_first(grid.T)
    grid = grid.reindex(index=ordered_aas, columns=ordered_aas)
    # Only keep bottom triangle
    grid = grid.where(np.tril(np.ones(grid.shape), k=0).astype(bool))
    plt.figure(figsize=(8, 6))
    cmap = plt.cm.viridis
    im = plt.imshow(grid, cmap=cmap, interpolation="nearest")
    plt.colorbar(im, label="Count")
    plt.xticks(range(len(grid.columns)), grid.columns)
    plt.yticks(range(len(grid.index)), grid.index)
    ax = plt.gca()
    sns.despine(ax=ax, top=True, right=True)
    plt.savefig(plot_path, dpi=600)


def find_specific_peptides(df, selected_peptide_list):
    found_peptides = (
        df[(df["similarity_score"] > 0.7)]["peptide_search"]
        # .str.replace("/1", "")
        .str.replace("/2", "").drop_duplicates()
    )
    found_isomers = (
        df[(df["similarity_score"] > 0.7)]["peptide_db"]
        .str.replace("/1", "")
        .drop_duplicates()
    )
    print(f"Found {len(found_peptides)} unique peptides with similarity > 0.7")
    print(f"Found {len(found_isomers)} unique isomers with similarity > 0.7")
    print(
        f"Selected peptides found in peptides: {len(found_peptides[found_peptides.isin(selected_peptide_list)])}"
    )
    print(
        f"Selected peptides found in isomers: {len(found_isomers[found_isomers.isin(selected_peptide_list)])}"
    )
    found_selected_peptides = []
    for pep in selected_peptide_list:
        if pep in found_peptides.values:
            found_selected_peptides.append(pep)
        else:
            continue
    return found_selected_peptides
