import pandas as pd
import numpy as np
from utils import read_pool_dir
from find_positional_isomers import find_positional_isomers, find_I_L_isomers
from positional_analysis import label_positional_isomers, plot_aa_grid_bottom_triangle

data_folder = "/Users/adams/Code/positional-isomer-distinguishability/data"

df = read_pool_dir(data_folder)
len(df["Sequence"].unique().tolist())
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
