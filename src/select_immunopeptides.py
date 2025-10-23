import os
import pandas as pd

path = "/Users/adams/Projects/indistinguishable-peptides/hla-query/SysteMHC-Atlas-non-UniProt.xlsx"
df = pd.read_excel(path, sheet_name="Blad1")

df.columns

# Add Peptide counts
peptide_counts = df["Peptide"].value_counts()
df["Peptide_Count"] = df["Peptide"].map(peptide_counts)

peptide_df = df[(df["Binder"].isin(["Yes", "No"])) & (df["Peptide_Count"] > 20)][
    ["Peptide"]
].drop_duplicates()

peptide_df["Peptide_Length"] = peptide_df["Peptide"].str.len()


peptide_df[
    (8 <= peptide_df["Peptide_Length"]) & (peptide_df["Peptide_Length"] <= 12)
].to_csv(
    "/Users/adams/Projects/indistinguishable-peptides/hla-query/systemhc_yesnobinder_20pepcount.csv",
    index=False,
    header=False,
)
