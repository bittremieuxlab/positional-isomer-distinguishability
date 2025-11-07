import pandas as pd
from matplotlib import pyplot as plt
import seaborn as sns
from matplotlib.ticker import FuncFormatter


def plot_db_proportions(plot_df: pd.DataFrame, plot_path: str = None):
    "Boxplot for all measurement at each diversity level"
    plt.figure(figsize=(9, 6))
    sns.set_context("talk")
    sns.set_style("whitegrid", rc={"grid.linewidth": 0.5})
    sns.boxplot(x="diversity", y="value", data=plot_df, color="#8ac6d0")
    sns.despine()
    plt.xlabel("Positional isomer fraction")
    plt.ylabel("Fraction with similarity score >= 0.7")
    plt.savefig(plot_path + "/db_proportions.png", dpi=300)


def plot_db_sizes(plot_df: pd.DataFrame, plot_path: str = None):
    "Boxplot for all measurement at each size"
    plt.figure(figsize=(9, 6))
    sns.set_context("talk")
    sns.set_style("whitegrid", rc={"grid.linewidth": 0.5})
    sns.boxplot(x="size", y="value", data=plot_df, color="#8ac6d0")
    plt.xlabel("Number of peptides")
    plt.ylabel("Fraction with similarity score >= 0.7")
    sizes = plot_df["size"].unique()
    sns.despine()
    plt.xticks(ticks=range(len(sizes)), labels=[f"{val//1000}K" for val in sizes])
    plt.savefig(plot_path + "/db_sizes.png", dpi=300)


def plot_db_sizes_proportions(plot_df: pd.DataFrame, plot_path: str = None):
    plt.figure(figsize=(10, 6))
    sns.set_context("talk")
    # colors = {"0.0": "#8ac6d0", "0.25": "#0e1c36", "0.5": "#2ca02c", "0.75": "#0a5239"}
    colors = {"0.0": "#8ac6d0"}
    for diversity, group in plot_df.groupby("diversity"):
        sns.lineplot(
            x=group["size"],
            y=group["value"],
            marker="o",
            label=f"Positional isomer fraction: {diversity}",
            color=colors[str(diversity)],
        )
    sns.despine()
    sns.set_style("whitegrid", rc={"grid.linewidth": 0.5})
    plt.xlabel("Number of peptides")
    plt.ylabel("Fraction with similarity score >= 0.7")

    # Format x-axis labels as 'K' for thousands
    def format_k(x, pos):
        if x >= 1000:
            return f"{int(x/1000)}K"
        return str(int(x))

    plt.gca().xaxis.set_major_formatter(FuncFormatter(format_k))
    plt.legend()
    plt.savefig(plot_path + "/db_sizes_proportions_0.png", dpi=300)


path = "/Users/adams/Projects/indistinguishable-peptides/db_size_diversity/db_sizes.csv"
plot_path = "/Users/adams/Projects/indistinguishable-peptides/Figures"
plot_df = pd.read_csv(path)

plot_db_proportions(plot_df, plot_path)
plot_db_sizes(plot_df, plot_path)
plot_db_sizes_proportions(plot_df, plot_path)
plot_db_sizes_proportions(plot_df[plot_df["diversity"] == 0], plot_path)
