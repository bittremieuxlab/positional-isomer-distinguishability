import pandas as pd
import numpy as np
import os
from generate_random_dbs import generate_random_db
from db_sizes_test import replace_with_positional_isomers
from find_positional_isomers import find_positional_isomers
from positional_analysis import label_positional_isomers
from koinapy import Koina
import requests

# df = generate_random_db(
#     dir_path="/Users/adams/Code/positional-isomer-distinguishability/data", size=None
# )
# div_df = replace_with_positional_isomers(df, frac=0.75)

inputs = pd.DataFrame()
# inputs["peptide_sequences"] = np.array(div_df["Sequence"])
inputs["peptide_sequences"] = [
    "KLEAGIR",
    "KIEAGLR",
    "RSNKPIILR",
    "RSNKPILIR",
]

div_df = pd.DataFrame()
div_df["Sequence"] = inputs["peptide_sequences"]

model = Koina("Prosit_2019_irt", "koina.wilhelmlab.org:443")
# model = Koina("Deeplc_hela_hf", "koina.wilhelmlab.org:443")
# model = Koina("Chronologer_RT", "koina.wilhelmlab.org:443")
predictions = model.predict(inputs)
div_df["predicted_irt"] = predictions["irt"].tolist()

isomers = find_positional_isomers(div_df["Sequence"].unique().tolist(), max_workers=8)
df_isomers = pd.DataFrame(isomers, columns=["Peptide1", "Peptide2"])

df_isomers = df_isomers.merge(
    div_df[["Sequence", "predicted_irt"]], left_on="Peptide1", right_on="Sequence"
).rename(columns={"predicted_irt": "Predicted_iRT_1"})
df_isomers = df_isomers.merge(
    div_df[["Sequence", "predicted_irt"]], left_on="Peptide2", right_on="Sequence"
).rename(columns={"predicted_irt": "Predicted_iRT_2"})
df_isomers = df_isomers.drop(columns=["Sequence_x", "Sequence_y"])
df_isomers["iRT_diff"] = (
    df_isomers["Predicted_iRT_1"] - df_isomers["Predicted_iRT_2"]
).abs()

df_isomers["peptide_length"] = df_isomers["Peptide1"].str.len()

df_isomers_labeled = label_positional_isomers(
    df_isomers, pep1="Peptide1", pep2="Peptide2"
)

# Plot iRT differences by peptide length
import matplotlib.pyplot as plt
import seaborn as sns

plt.figure(figsize=(10, 6))
sns.boxplot(
    data=df_isomers_labeled,
    x="peptide_length",
    y="iRT_diff",
    showfliers=False,
)
plt.title("iRT Differences of Positional Isomers by Peptide Length")
plt.xlabel("Peptide Length")
plt.ylabel("Absolute iRT Difference")
plt.grid(True, which="both", linestyle="--", linewidth=0.5)
plt.show()

# Add value counts to df_isomers_labeled
df_isomers_labeled["positional_aa_counts"] = df_isomers_labeled["positional_aa"].map(
    df_isomers_labeled["positional_aa"].value_counts()
)

df_isomers_labeled = df_isomers_labeled.sort_values(by="positional_aa")
plt.figure(figsize=(10, 6))
sns.boxplot(
    data=df_isomers_labeled[df_isomers_labeled["positional_aa_counts"] < 1000],
    x="positional_aa",
    y="iRT_diff",
    showfliers=False,
)
plt.title("Similarity Score by Positional Isomer Amino Acids")
plt.xlabel("Positional Isomer Amino Acids")
plt.ylabel("Absolute iRT Difference")
plt.xticks(rotation=45)
plt.grid(axis="y")
plt.show()

plt.close()
