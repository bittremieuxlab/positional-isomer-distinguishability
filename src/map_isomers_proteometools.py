import pandas as pd
from src.utils import read_pool_dir
from src.proteome_tools_peptides import find_positional_isomers


data_folder = "data"

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

# Map Pool_name back to df_isomers
df_merged = df_isomers.merge(
    df[["Sequence", "Pool_name"]], left_on="Peptide1", right_on="Sequence"
).merge(
    df[["Sequence", "Pool_name"]],
    left_on="Peptide2",
    right_on="Sequence",
    suffixes=("_1", "_2"),
)
df_merged = df_merged[["Peptide1", "Peptide2", "Pool_name_1", "Pool_name_2"]]

df_merged.value_counts("Pool_name_1")
df_merged.value_counts("Pool_name_2")


df_merged[df_merged["Pool_name_1"] == df_merged["Pool_name_2"]].value_counts(
    "Pool_name_1"
)
df_merged[df_merged["Pool_name_1"].str.contains("HLA")].value_counts("Pool_name_1")

df_merged["lenth"] = df_merged["Peptide1"].apply(len)
df_merged.value_counts("lenth")

df_merged["last_aa_1"] = df_merged["Peptide1"].str[-1]
df_merged["last_aa_2"] = df_merged["Peptide2"].str[-1]

df_merged[
    df_merged["Pool_name_1"].str.contains("_HLA_")
    & df_merged["Pool_name_2"].str.contains("_HLA_")
].value_counts("last_aa_1")

df_merged[
    df_merged["Pool_name_1"].str.contains("_HLA_")
    & df_merged["Pool_name_2"].str.contains("_HLA_")
].value_counts("last_aa_2")


df_merged.columns

len(df_merged)

df_merged.to_csv(
    "positional_isomers_msci_input_mapped.csv",
    index=False,
)
