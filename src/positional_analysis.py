import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
import scipy.stats as stats


def is_positional_isomer(peptide1: str, peptide2: str) -> bool:
    return peptide1 != peptide2 and sorted(peptide1) == sorted(peptide2)


def label_positional_isomers(df):
    is_positional_isomer_list = []
    positional_isomer_distance_list = []
    positional_aa_list = []
    first_difference_position_list = []
    second_difference_position_list = []
    length_list = []
    for peptide1, peptide2 in zip(df["peptide 1"], df["peptide 2"]):
        diffs = [(a, b) for a, b in zip(peptide1, peptide2) if a != b]
        is_positional_isomer_list.append(len(diffs) == 2 and diffs[0] == diffs[1][::-1])
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
    return df


def plot_isomer_distance(df):
    plt.figure(figsize=(10, 6))
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
    plt.show()


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
    plt.figure(figsize=(10, 6))
    sns.boxplot(
        data=df,
        x="positional_aa",
        y="similarity_score",
    )
    plt.title("Similarity Score by Positional Isomer Amino Acids")
    plt.xlabel("Positional Isomer Amino Acids")
    plt.ylabel("Similarity Score")
    plt.ylim(0, 1)
    plt.xticks(rotation=45)
    plt.grid(axis="y")
    plt.show()


msci_output_path = "/Users/adams/Projects/indistinguishable-peptides/positional-isomers/positional_isomers_msci_output.csv"

df = pd.read_csv(msci_output_path)
df = label_positional_isomers(df)

df["RT_diff"] = np.abs(df["iRT 1"] - df["iRT 2"])
df["positional_amino_acids"] = df["positional_amino_acids"].astype(str)

# concat('I', 'V') to IV,  ('S', 'G') to SG, etc. and sort aa alphabetically
df["positional_aa"] = (
    df["positional_amino_acids"]
    .apply(lambda x: "".join(eval(x)))
    .apply(lambda x: "".join(sorted(x)))
)
# Order df based on positional_aa
df = df.sort_values(by="positional_aa")

plot_isomer_position(
    df[
        (df["is_positional_isomer"] == True)
        & (df["is_positional_isomer_distance"] == 1)
        & (df["RT_diff"] < 0.5)
    ]
)

plot_isomer_aa(
    df[
        # (df["similarity_score"] > 0.7)
        (df["is_positional_isomer"] == True)
        # & (df["is_positional_isomer_distance"] == 1)
        # (df["RT_diff"] < 0.5)
        & (df["RT_diff"] < 0.5)
    ]
)
