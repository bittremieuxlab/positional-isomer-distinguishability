import pandas as pd
import numpy as np
from utils import read_pool_dir
from find_positional_isomers import find_positional_isomers, find_I_L_isomers
from positional_analysis import (
    label_positional_isomers,
    plot_aa_grid_bottom_triangle,
    plot_position_length_ratio_vs_similarity,
    plot_AA_length_grid_bottom_triangle,
)
import matplotlib.pyplot as plt
import seaborn as sns

data_folder = "/Users/adams/Code/positional-isomer-distinguishability/data"

df = read_pool_dir(data_folder)
len(df["Sequence"].unique().tolist())

# Look for peptides with "YW" or "WY" in the sequence
df_yw = df[df["Sequence"].str.contains("YW|WY")]
len(df_yw["Sequence"].unique().tolist())

# Generate a list of positional isomers by swapping Y and W
isomers_yw = []
for seq in df_yw["Sequence"].unique().tolist():
    if "YW" in seq:
        isomer_seq = seq.replace("YW", "WY")
        isomers_yw.append((seq, isomer_seq))
    elif "WY" in seq:
        isomer_seq = seq.replace("WY", "YW")
        isomers_yw.append((seq, isomer_seq))

# Create an input for MSCI (so single column with sequence)
file_path = "/Users/adams/Projects/indistinguishable-peptides/proteometools_yw_pep1.csv"
with open(file_path, "w") as f:
    for seq1, seq2 in isomers_yw:
        f.write(f"{seq1}\n")

file_path = "/Users/adams/Projects/indistinguishable-peptides/proteometools_yw_pep2.csv"
with open(file_path, "w") as f:
    for seq1, seq2 in isomers_yw:
        f.write(f"{seq2}\n")

msci_output_path = "/Users/adams/Projects/indistinguishable-peptides/proteometools_yw_db_charge_2_deeplc.csv"
df_yw_msci = pd.read_csv(msci_output_path)
df_yw_msci = label_positional_isomers(
    df_yw_msci, pep1="peptide_search", pep2="peptide_db"
)
df_yw_msci.columns
df_yw_msci[(df_yw_msci["positional_aa"] == "WY")][
    ["peptide_search", "peptide_db", "similarity_score"]
]

len(df_yw_msci[(df_yw_msci["positional_aa"] == "WY")])
df_yw_msci[(df_yw_msci["positional_aa"] == "WY")]["similarity_score"].describe()
df_yw_msci["peptide_search"] = df_yw_msci["peptide_search"].str.replace("/2", "")
df_yw_msci["peptide_db"] = df_yw_msci["peptide_db"].str.replace("/2", "")
df_yw_msci["peptide_length"] = df_yw_msci["peptide_search"].apply(len)
df_yw_msci["r_count"] = df_yw_msci["peptide_search"].str.count("R")
df_yw_msci["p_count"] = df_yw_msci["peptide_search"].str.count("P")

df_yw_msci["position_length_ratio_start"] = (
    df_yw_msci["first_difference_position"] / df_yw_msci["peptide_length"]
)
df_yw_msci["position_length_ratio_end"] = (
    df_yw_msci["second_difference_position"] / df_yw_msci["peptide_length"]
)

for r in range(max(df_yw_msci["r_count"]) + 1):
    # plot_isomer_position(
    #     df_yw_msci[
    #         (df_yw_msci["positional_aa"] == "WY")
    #         # & (df_yw_msci["peptide_length"] == 9)
    #         & (df_yw_msci["is_positional_isomer_distance"] == 1)
    #         & (df_yw_msci["r_count"] == r)
    #     ],
    #     plot_title=f"Proteometools YW isomers with {r} R residues",
    #     # f"/Users/adams/Projects/indistinguishable-peptides/Figures/proteometools_yw_isomer_position_9_{r}R.png",
    #     plot_path=f"/Users/adams/Projects/indistinguishable-peptides/Figures/proteometools_yw_isomer_{r}R.png",
    # )
    for position_length_ratio in [
        "position_length_ratio_start",
        "position_length_ratio_end",
    ]:
        plot_position_length_ratio_vs_similarity(
            df_yw_msci[
                (df_yw_msci["positional_aa"] == "WY")
                & (df_yw_msci["is_positional_isomer_distance"] == 1)
                & (df_yw_msci["r_count"] == r)
            ],
            position_length_ratio=position_length_ratio,
            plot_path=f"/Users/adams/Projects/indistinguishable-peptides/Figures/proteometools_yw{position_length_ratio}_vs_similarity_{r}R.png",
            plot_title=f"Proteometools YW Isomers - {position_length_ratio.replace('_', ' ').title()} - {r} R residues",
        )


plot_AA_length_grid_bottom_triangle(
    df_yw_msci[
        (df_yw_msci["positional_aa"] == "WY")
        & (df_yw_msci["is_positional_isomer_distance"] == 1)
    ],
    count_column="p_count",
    count_name="P Count",
    plot_path=f"/Users/adams/Projects/indistinguishable-peptides/Figures/proteometools_yw_P_count_grid.png",
)


df_yw_msci[
    (df_yw_msci["positional_aa"] == "WY")
    & (df_yw_msci["is_positional_isomer_distance"] == 1)
    & (df_yw_msci["peptide_length"] == 9)
][["peptide_search", "peptide_db", "similarity_score", "first_difference_position"]]


df_yw_msci[
    (df_yw_msci["positional_aa"] == "WY")
    # & (df_yw_msci["peptide_length"] == 9)
    & (df_yw_msci["is_positional_isomer_distance"] == 1)
    & (df_yw_msci["r_count"] == 3)
][["peptide_search", "peptide_db", "similarity_score", "position_length_ratio"]]

df_yw_msci[(df_yw_msci["positional_aa"] == "WY")].peptide_length.value_counts()

# Investigate positional isomers in proteometools data
isomers = find_positional_isomers(df["Sequence"].unique().tolist(), max_workers=8)
df_isomers = pd.DataFrame(isomers, columns=["Peptide1", "Peptide2"])
df_isomers = label_positional_isomers(df_isomers, pep1="Peptide1", pep2="Peptide2")

IL_isomers = find_I_L_isomers(df["Sequence"].unique().tolist())
df_IL_isomers = pd.DataFrame(IL_isomers, columns=["Peptide1", "Peptide2"])

overlapping_isomers = df_isomers[
    (df_isomers["Peptide1"] + df_isomers["Peptide2"]).isin(
        df_IL_isomers["Peptide1"] + df_IL_isomers["Peptide2"]
    )
]

len(overlapping_isomers)
len(df_isomers)

df_isomers.columns
len(df_isomers[df_isomers["is_positional_isomer_distance"] == 1])


# plot AA grid showing which AAs are involved in positional isomers in the grid show the count

plot_path = "/Users/adams/Projects/indistinguishable-peptides/Figures"
plot_aa_grid_bottom_triangle(
    df_isomers,
    # df_isomers[df_isomers["is_positional_isomer_distance"] == 1],
    f"{plot_path}/positional_isomer_aa_grid.png",
)

df_isomers[["aa1", "aa2"]] = df_isomers["positional_aa"].apply(
    lambda x: pd.Series(list(x))
)

counts = df_isomers.groupby(["aa1", "aa2"]).size().reset_index(name="count")
df_isomers.value_counts("positional_aa").nlargest(20)


excel_path = "/Users/adams/Downloads/41467_2021_23713_MOESM8_ESM.xlsx"

# ignore first line
excel_df = pd.read_excel(excel_path, sheet_name="Combined results", skiprows=1)
excel_df.columns

# Order the aa in the peptide sequence alphabetically to identify positional isomers
excel_df["peptide_sequence_content"] = excel_df["Peptide Sequence"].apply(
    lambda x: "".join(sorted(list(x)))
)
excel_df["sequence_content"] = (
    excel_df["Sequence"]
    .astype(str)
    .apply(lambda x: "".join(sorted(x)) if x != "nan" else None)
)

excel_df[
    (excel_df["Peptide type"] == "spliced")
    & (excel_df["peptide_sequence_content"] == excel_df["sequence_content"])
    & (excel_df["Spectral angle"] > 0.7)
    & (excel_df["No better hypothesis with MSFragger"] == "yes")
]
