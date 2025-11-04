import pandas as pd
from positional_analysis import (
    label_positional_isomers,
    plot_isomer_position,
    plot_isomer_aa,
)
from utils import report_isomer_stats
import numpy as np


msci_output_path1 = "/Users/adams/Projects/indistinguishable-peptides/positional-isomers/positional_isomers_msci_output.csv"
msci_output_path2 = "/Users/adams/Projects/indistinguishable-peptides/positional-isomers/positional_isomers_msci_output_tims.csv"
# msci_output_path = "/Users/adams/Projects/indistinguishable-peptides/hla-query/hla_noPTM_peptides_canonical_human_peptides_8_12_output.csv"
# msci_output_path = "/Users/adams/Projects/indistinguishable-peptides/hla-query/hla_noPTM_peptides_canonical_human_peptides_8_12_charge_1.csv"
msci_output_path = "/Users/adams/Projects/indistinguishable-peptides/hla-query/charlottes13peptides_canonical_human_peptides_8_12_charge_2.csv"
msci_output_path = "/Users/adams/Projects/indistinguishable-peptides/hla-query/mhcquant3peptides_canonical_human_peptides_8_12_charge_2.csv"
msci_output_path = "/Users/adams/Projects/indistinguishable-peptides/hla-query/systemhc_yesnobinder_20pepcount_charge_1_rt10.csv"

df = pd.read_csv(msci_output_path)
df = label_positional_isomers(df, pep1="peptide_search", pep2="peptide_db")

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

df["rt_diff"] = np.abs(df["iRT_search"] - df["iRT_db"])
df[(df["similarity_score"] > 0.7)].rt_diff

df.similarity_score.max()

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

report_isomer_stats(df)
plot_isomer_position(
    df[
        (df["is_positional_isomer"] == True)
        & (df["is_positional_isomer_distance"] == 1)
        # & (df["RT_diff"] < 0.5)
    ]
)
plot_isomer_aa(
    df[
        (df["similarity_score"] > 0.7)
        & (df["peptide_search"] != df["peptide_db"])
        # (df["is_positional_isomer"] == True)
        # & (df["is_positional_isomer_distance"] == 1)
        # (df["RT_diff"] < 0.5)
        # & (df["RT_diff"] < 0.5)
    ]
)
found_peptides_df = df[(df["similarity_score"] > 0.7)]
found_peptides_df["found_peptide_search"] = found_peptides_df[
    "peptide_search"
].str.replace("/2", "")
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
found_peptides_df["found_peptide_search"].isin(selected_peptide_list)
found_peptides_df[
    found_peptides_df["peptide_search"] != found_peptides_df["peptide_db"]
]
