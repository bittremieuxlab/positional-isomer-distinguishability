import pandas as pd
from positional_analysis import (
    label_positional_isomers,
    plot_AA_length_grid_bottom_triangle,
    # plot_isomer_distance,
    # plot_isomer_aa,
)
import seaborn as sns
from matplotlib.ticker import FuncFormatter
from matplotlib import pyplot as plt
import numpy as np
from koinapy import Koina
from scipy.stats import mannwhitneyu


def plot_isomer_distance(df, plot_path: str = None):
    plt.figure(figsize=(12, 6))
    sns.set_context("talk")
    ax = sns.boxplot(
        data=df,
        x="is_positional_isomer_distance",
        y="similarity_score",
    )
    # plt.title("Similarity Score by Positional Isomer Distance")
    plt.xlabel("Positional Isomer Distance")
    plt.ylabel("Similarity Score")
    plt.ylim(0, 1)
    plt.grid(axis="y")
    sns.despine()
    unique_distances = sorted(df["is_positional_isomer_distance"].unique())
    ax.set_xticks(range(len(unique_distances)))
    ax.set_xticklabels([str(int(x)) for x in unique_distances])
    # plt.tight_layout()
    plt.savefig(plot_path + "/isomer_distance.png", dpi=300)


def plot_isomer_distance_rt(df, plot_path: str = None):
    plt.figure(figsize=(12, 6))
    sns.set_context("talk")
    ax = sns.boxplot(
        data=df,
        x="is_positional_isomer_distance",
        y="rt_diff",
    )
    # plt.title("iRT Difference by Positional Isomer Distance")
    plt.xlabel("Positional Isomer Distance")
    plt.ylabel("iRT Difference")
    plt.grid(axis="y")
    sns.despine()
    unique_distances = sorted(df["is_positional_isomer_distance"].unique())
    ax.set_xticks(range(len(unique_distances)))
    ax.set_xticklabels([str(int(x)) for x in unique_distances])
    # plt.tight_layout()
    plt.savefig(plot_path, dpi=300)


def plot_aa_grid(
    df,
    median_value="similarity_score",
    median_name="Median Similarity",
    plot_path: str = None,
):
    # Define amino acid groups by property
    aa_groups = {
        "Nonpolar": ["A", "V", "L", "I", "M", "F", "W", "P", "G"],
        "Polar (uncharged)": ["S", "T", "C", "Y", "N", "Q"],
        "Positive (basic)": ["K", "R", "H"],
        "Negative (acidic)": ["D", "E"],
    }
    ordered_aas = [aa for group in aa_groups.values() for aa in group]
    # Split into separate AAs
    df[["aa1", "aa2"]] = df["positional_aa"].apply(lambda x: pd.Series(list(x)))
    # Median similarity per unique pair
    median_sim = df.groupby(["aa1", "aa2"])[median_value].median().reset_index()
    # Pivot into symmetric matrix
    grid = median_sim.pivot(index="aa1", columns="aa2", values=median_value)
    # Ensure full symmetric matrix (fill both triangles)
    grid = grid.combine_first(grid.T)
    grid = grid.reindex(index=ordered_aas, columns=ordered_aas)
    plt.figure(figsize=(8, 6))
    cmap = plt.cm.viridis
    im = plt.imshow(grid, cmap=cmap, interpolation="nearest")
    cbar = plt.colorbar(im, label=median_name)
    # Annotate only visible values
    # for (i, j), val in np.ndenumerate(grid.values):
    #     plt.text(j, i, f"{val:.2f}", ha="center", va="center", color="w")
    plt.grid(visible=False)
    plt.tick_params(bottom=True, left=True)
    plt.xticks(range(len(grid.columns)), grid.columns)
    plt.yticks(range(len(grid.index)), grid.index)
    # Get the current Axes
    ax = plt.gca()
    # Remove internal gridlines and make outer border black
    plt.grid(False)
    for spine in ax.spines.values():
        spine.set_visible(True)
        spine.set_color("black")
        spine.set_linewidth(1.5)
    for spine in cbar.ax.spines.values():
        spine.set_visible(True)
        spine.set_color("black")
        spine.set_linewidth(1.5)
    plt.tick_params(bottom=True, left=True, color="black", labelcolor="black")
    plt.savefig(plot_path, dpi=300)


def plot_peptide_length(
    df,
    plot_path: str = None,
):
    plt.figure(figsize=(12, 6))
    sns.set_context("talk")
    df["peptide_length"] = df["peptide_length"].astype(int)
    sns.histplot(data=df, x="peptide_length", discrete=True, color="#8ac6d0")
    plt.xlim(6, 31)
    plt.xlabel("Peptide length")
    plt.ylabel("Amount")
    plt.grid(axis="y")
    sns.despine()
    plt.savefig(plot_path, dpi=300)


def plot_rt_pep1_pep2(df, plot_path: str = None, min=-50, max=199):
    plt.figure(figsize=(8, 8))
    sns.set_context("talk")
    sns.scatterplot(data=df, x="iRT 1", y="iRT 2", alpha=0.5, color="#8ac6d0")
    plt.xlabel("iRT Peptide 1")
    plt.ylabel("iRT Peptide 2")
    plt.xlim(min, max)
    plt.ylim(min, max)
    sns.despine()
    plt.savefig(plot_path, dpi=300)


msci_output_path = (
    # "/Users/adams/Projects/indistinguishable-peptides/db_size_diversity/300000_0.25.csv"
    "/Users/adams/Projects/indistinguishable-peptides/db_size_diversity/100000_0.5_50.csv"
    # "/Users/adams/Projects/indistinguishable-peptides/positional-isomers/100000_0.5.csv"
)

df = pd.read_csv(msci_output_path)
df = label_positional_isomers(df, pep1="peptide 1", pep2="peptide 2")

df["rt_diff"] = (df["iRT 1"] - df["iRT 2"]).abs()
df["peptide 1"] = df["peptide 1"].str.replace("/2", "")
df["peptide 2"] = df["peptide 2"].str.replace("/2", "")
df["peptide_length"] = df["peptide 1"].apply(len)
df["position_length_ratio"] = df["first_difference_position"] / df["peptide_length"]
df["position_length_ratio_start"] = (
    df["first_difference_position"] / df["peptide_length"]
)
df["position_length_ratio_end"] = (
    df["second_difference_position"] / df["peptide_length"]
)
df["R_count"] = df["peptide 1"].str.count("R")
df["P_count"] = df["peptide 1"].str.count("P")

df[
    (df["peptide 1"] != df["peptide 2"])
    & (df["is_positional_isomer"] == True)
    # & (df["is_positional_isomer_distance"] == 1)
    # & (df["P_count"] == 0)
    # & (df["positional_aa"] == "IL")
    # & (df["positional_aa"].str.contains("P"))
    & (df["positional_aa"].str.contains("P") == False)
]["similarity_score"].median()

df[
    (df["peptide 1"] != df["peptide 2"])
    & (df["is_positional_isomer"] == True)
    # & (df["is_positional_isomer_distance"] == 1)
    & (df["position_length_ratio_start"] < 0.2)
    # & (df["position_length_ratio_end"] > 0.8)
    # & (df["position_length_ratio"] > 0.8)
]["similarity_score"].median()


# plot_path_grid = "/Users/adams/Projects/indistinguishable-peptides/Figures/positional_isomer_R_length_grid.png"
plot_path_grid = "/Users/adams/Projects/indistinguishable-peptides/Figures/positional_isomer_P_length_grid.png"
plot_AA_length_grid_bottom_triangle(
    df[
        (df["peptide 1"] != df["peptide 2"])
        & (df["is_positional_isomer"] == True)
        & (df["is_positional_isomer_distance"] == 1)
    ],
    count_column="P_count",
    count_name="P Count",
    plot_path=plot_path_grid,
)


plot_path = "/Users/adams/Projects/indistinguishable-peptides/Figures"

plot_isomer_distance(
    df[
        (df["is_positional_isomer"] == True)
        # & (df["is_positional_isomer_distance"] == 1)
        # & (df["RT_diff"] < 0.5)
    ],
    plot_path,
)

plot_isomer_distance_rt(
    df[
        (df["is_positional_isomer"] == True)
        # & (df["similarity_score"] > 0.7)
    ],
    # plot_path=plot_path + "/isomer_distance_rt_simcutoff.png",
    plot_path=plot_path + "/isomer_distance_rt.png",
)

plot_peptide_length(
    df[
        (df["is_positional_isomer"] == True)
        & (df["similarity_score"] > 0.7)
        # & (df["is_positional_isomer_distance"] < 5)
        & (df["rt_diff"] > 20)
    ],
    # plot_path=plot_path + "/peptide_length_sim07.png",
    plot_path=plot_path + "/peptide_length_sim07_rt20.png",
)

plot_rt_pep1_pep2(
    df_chronologer[
        (df_chronologer["is_positional_isomer"] == True)
        & (df_chronologer["positional_aa"] != "IL")
        & (df_chronologer["is_positional_isomer_distance"] == 1)
    ],
    plot_path=plot_path + "/rt_chronologer_pep1_pep2_notIL_dist1.png",
    min=-2,
    max=29,
    # plot_path=plot_path + "/rt_pep1_pep2_notIL.png",
)

inputs = pd.DataFrame()
inputs["peptide_sequences"] = np.array(
    df[df["is_positional_isomer"] == True]["peptide 1"].str.replace("/2", "")
)
model = Koina("Chronologer_RT", "koina.wilhelmlab.org:443")
predictions = model.predict(inputs)
# df_chronologer = df[df["is_positional_isomer"] == True]
df_chronologer["iRT 1"] = predictions["rt"].tolist()

df_chronologer["rt_diff"] = (df_chronologer["iRT 1"] - df_chronologer["iRT 2"]).abs()
df_deeplc["rt_diff"] = (df_deeplc["iRT 1"] - df_deeplc["iRT 2"]).abs()
df["rt_diff"] = (df["iRT 1"] - df["iRT 2"]).abs()

df_chronologer[
    (df_chronologer["is_positional_isomer"] == True)
    & (df_chronologer["positional_aa"] == "IK")
    & (df_chronologer["similarity_score"] > 0.7)
][["peptide 1", "peptide 2", "iRT 1", "iRT 2", "rt_diff"]]

df_deeplc[
    (df_deeplc["is_positional_isomer"] == True)
    & (df_deeplc["positional_aa"] == "IK")
    & (df_deeplc["similarity_score"] > 0.7)
][["peptide 1", "peptide 2", "iRT 1", "iRT 2", "rt_diff"]]

plot_isomer_distance_rt(
    df_deeplc[
        (df_deeplc["is_positional_isomer"] == True)
        & (df_deeplc["similarity_score"] > 0.7)
    ],
    plot_path=plot_path + "/isomer_distance_rt_simcutoff_deeplc.png",
)

df_deeplc[
    (df_deeplc["is_positional_isomer"] == True) & (df_deeplc["similarity_score"] > 0.7)
]["rt_diff"].median()

df[
    (df["is_positional_isomer"] == True)
    # & (df["positional_aa"] != "IL")
    # & (df["rt_diff"] < 5)
    # & (df["is_positional_isomer_distance"] == 1)
    & (df["similarity_score"] > 0.7)
    & (df["positional_aa"] == "IK")
][
    [
        "peptide 1",
        "peptide 2",
        "iRT 1",
        "iRT 2",
        "rt_diff",
        "is_positional_isomer_distance",
    ]
]

df[
    (df["is_positional_isomer"] == True)
    # & (df["similarity_score"] > 0.7)
    & (df["rt_diff"] > 30)
    & (df["is_positional_isomer_distance"] < 5)
    # & (df["is_positional_isomer_distance"] == 9)
]["peptide_length"].describe()

df[
    (df["is_positional_isomer"] == True)
    & (df["similarity_score"] > 0.7)
    & (df["rt_diff"] > 20)
    & (df["peptide_length"] > 7)
    & (df["is_positional_isomer_distance"] > 5)
][["peptide 1", "peptide 2", "positional_aa", "peptide_length"]].value_counts(
    "positional_aa"
)


plot_aa_grid(
    df[
        (df["is_positional_isomer"] == True)
        # & (df["similarity_score"] > 0.7)
        & (df["is_positional_isomer_distance"] == 1)
    ],
    median_value="rt_diff",
    median_name="Median iRT Difference",
    plot_path=plot_path + "/isomer_aa_grid_rt_dist1.png",
)

# median rt difference
df[
    (df["is_positional_isomer"] == True)
    # & (df["similarity_score"] > 0.7)
]["rt_diff"].mean()


df[(df["is_positional_isomer"] == True)]["is_positional_isomer_distance"].value_counts()


# Add value counts to df_isomers_labeled
df["positional_aa_counts"] = df["positional_aa"].map(df["positional_aa"].value_counts())
