import pandas as pd

query_path = (
    "/Users/adams/Projects/indistinguishable-peptides/hla-query/hla_noPTM_peptides.csv"
)
db_path = "/Users/adams/Projects/indistinguishable-peptides/hla-query/canonical_human_peptides_8_12.csv"
query_df = pd.read_csv(query_path, header=None, names=["Peptide"])
db_df = pd.read_csv(db_path, header=None, names=["Peptide"])

db_df["IL_peptide"] = db_df["Peptide"].str.replace("I", "L")
query_df["IL_peptide"] = query_df["Peptide"].str.replace("I", "L")

db_df["aa_content"] = db_df["Peptide"].apply(lambda x: "".join(sorted(x)))
query_df["aa_content"] = query_df["Peptide"].apply(lambda x: "".join(sorted(x)))

db_df["IL_aa_content"] = db_df["IL_peptide"].apply(lambda x: "".join(sorted(x)))
query_df["IL_aa_content"] = query_df["IL_peptide"].apply(lambda x: "".join(sorted(x)))

len(set(db_df["IL_aa_content"]))
len(set(query_df["IL_aa_content"]))
len(set(db_df["aa_content"]))
len(set(query_df["aa_content"]))

# Filter query_df to only those with aa_content in db_df
filtered_db_df = (
    db_df[db_df["IL_aa_content"].isin(query_df["IL_aa_content"].unique())]
    .drop_duplicates()
    .reset_index(drop=True)
)

filtered_query_df = (
    query_df[query_df["IL_aa_content"].isin(db_df["IL_aa_content"].unique())]
    .drop_duplicates()
    .reset_index(drop=True)
)


len(filtered_db_df)
filtered_db_df["Peptide"].to_csv(
    "/Users/adams/Projects/indistinguishable-peptides/hla-query/filtered_canonical_human_peptides_8_12.csv",
    index=False,
    header=False,
)

path = "/Users/adams/Projects/indistinguishable-peptides/hla-query/filtered_canonical_human_peptides_8_12.csv"
filtered_db_df = pd.read_csv(path, header=None, names=["Peptide"])


# How many of the peptide sequences in filtered_db_df are also in query_df
len(set(query_df["Peptide"]).intersection(set(filtered_db_df["Peptide"])))
filtered_query_df = query_df[
    ~(query_df["IL_aa_content"].isin(filtered_db_df["IL_aa_content"].unique()))
]

filtered_query_df["length"] = filtered_query_df["Peptide"].str.len()
filtered_query_df["length"].value_counts()

len(set(query_df["IL_aa_content"]).intersection(set(filtered_db_df["IL_aa_content"])))

len(set(filtered_query_df["Peptide"]).intersection(set(filtered_db_df["Peptide"])))


selected_peptide_list = [
    "HELEKKIEL",
    "TFKEKVTSL",
    "VEEDFQQKL",
    "YEIDDVERL",
    "LPLHKVKSL",
    "KVAELLEKY",
    "DRYLLVSQF",
    "VVHLIKNAY",
    "TYSEKTTLF",
    "ALKARTVTF",
    "KLFTSTGLK",
    "KYLTIYLQK",
    "TPHPSELKVM",
]

spliced_peptides = set(selected_peptide_list).difference(set(filtered_db_df["Peptide"]))

spliced_peptide_df = pd.DataFrame(spliced_peptides, columns=["Peptide"])
spliced_peptide_df["IL_peptide"] = spliced_peptide_df["Peptide"].str.replace("I", "L")
spliced_peptide_df["IL_aa_content"] = spliced_peptide_df["IL_peptide"].apply(
    lambda x: "".join(sorted(x))
)

filtered_db_df["IL_peptide"] = filtered_db_df["Peptide"].str.replace("I", "L")
filtered_db_df["IL_aa_content"] = filtered_db_df["IL_peptide"].apply(
    lambda x: "".join(sorted(x))
)

spliced_db_peptides_df = filtered_db_df[
    filtered_db_df["IL_aa_content"].isin(spliced_peptide_df["IL_aa_content"])
]

spliced_df = spliced_peptide_df.merge(
    spliced_db_peptides_df, on="IL_aa_content", how="left", suffixes=("_spliced", "_db")
)[["Peptide_spliced", "Peptide_db"]]

spliced_df["diffs"] = [
    (a, b)
    for a, b in zip(spliced_df["Peptide_spliced"], spliced_df["Peptide_db"])
    if a != b
]


all_peptides_df = query_df.merge(
    filtered_db_df, on="IL_aa_content", how="left", suffixes=("_query", "_db")
)[["Peptide_query", "Peptide_db"]]

all_peptides_df.isna().sum()
# Remove rows where Peptide_db is NaN
all_peptides_df = all_peptides_df[~all_peptides_df["Peptide_db"].isna()]
# Remove where Peptide_query == Peptide_db
all_peptides_df = all_peptides_df[
    all_peptides_df["Peptide_query"] != all_peptides_df["Peptide_db"]
]
is_positional_isomer_list = []
positional_aa_list = []
for peptide1, peptide2 in zip(
    all_peptides_df["Peptide_query"], all_peptides_df["Peptide_db"]
):
    diffs = [(a, b) for a, b in zip(peptide1, peptide2) if a != b]
    is_positional_isomer_list.append(len(diffs) == 2 and diffs[0] == diffs[1][::-1])
    positional_aa_list.append(diffs[0])


is_positional_isomer_list.count(True)
all_peptides_df["is_positional_isomer"] = is_positional_isomer_list
all_peptides_df["positional_amino_acids"] = positional_aa_list

all_peptides_df[all_peptides_df["is_positional_isomer"]].value_counts(
    "positional_amino_acids"
)

all_peptides_df["len_diffs"] = all_peptides_df["diffs"].apply(
    lambda x: (len(x[0]), len(x[1]))
)
