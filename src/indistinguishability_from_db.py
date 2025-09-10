import numpy as np
import pandas as pd
from MSCI.Preprocessing.Parsing import read_msp_file

from msci_runner import (
    parallel_process_search_db_pairs,
)
from utils import (
    predict_spectra,
    process_peptide_db_combinations,
)
from pyteomics import mass

from pandarallel import pandarallel

pandarallel.initialize(progress_bar=True)


def filter_mz(
    search_pep_file: str,
    db_file: str,
    charge: int,
    mz_tolerance: int | float,
    out_db_file: str = None,
):
    search_df = pd.read_csv(search_pep_file, header=None)
    db_df = pd.read_csv(db_file, header=None)
    db_df = db_df[~db_df[0].str.contains("X")].reset_index(drop=True)

    search_df["m/z"] = search_df[0].parallel_apply(
        lambda x: mass.calculate_mass(sequence=x) / charge
    )
    db_df["m/z"] = db_df[0].parallel_apply(
        lambda x: mass.calculate_mass(sequence=x) / charge
    )

    mask = db_df["m/z"].parallel_apply(
        lambda mz: np.any(np.isclose(mz, search_df["m/z"], atol=mz_tolerance))
    )

    filtered_db_df = db_df[mask].reset_index(drop=True)
    if out_db_file is None:
        out_db_file = f"{db_file}_mz_filtered.csv"

    filtered_db_df.to_csv(out_db_file, columns=[0], index=False, header=False)
    print(f"Filtered db has {len(filtered_db_df)} peptides")

    return out_db_file


def indistinguishability_from_db(
    search_pep_file,
    db_file,
    collision_energy: int = 30,
    charge: int = 2,
    model_intensity: str = "Prosit_2020_intensity_HCD",
    model_irt: str = "Prosit_2019_irt",
    mz_tolerance: int | float = 1,
    irt_tolerance: int | float = 5,
    peak_tolerance: int | float = 0,
    peak_ppm: int | float = 10,
    output_file: str = "",
    n_chunks: int = None,
    prefilter_mz: bool = True,
):
    if prefilter_mz:
        db_file = filter_mz(search_pep_file, db_file, charge, mz_tolerance)

    search_spectra, search_pred_file = predict_spectra(
        input_file=search_pep_file,
        collision_energy=collision_energy,
        charge=charge,
        model_intensity=model_intensity,
        model_irt=model_irt,
    )

    db_spectra, db_pred_file = predict_spectra(
        input_file=db_file,
        collision_energy=collision_energy,
        charge=charge,
        model_intensity=model_intensity,
        model_irt=model_irt,
    )

    mz_irt_search_df = read_msp_file(search_pred_file)
    mz_irt_db_df = read_msp_file(db_pred_file)

    groups_df = process_peptide_db_combinations(
        mz_irt_search_df, mz_irt_db_df, mz_tolerance, irt_tolerance, use_ppm=False
    )

    index_array = groups_df[["index_search", "index_db"]].values.astype(int)
    result = parallel_process_search_db_pairs(
        index_array,
        n_chunks=n_chunks,
        spectra_search=search_spectra,
        spectra_db=db_spectra,
        mz_irt_search_df=mz_irt_search_df,
        mz_irt_db_df=mz_irt_db_df,
        tolerance=peak_tolerance,
        ppm=peak_ppm,
    )

    result.to_csv(output_file, index=False)


if __name__ == "__main__":
    search_pep_file = "search_db/charlottes13peptides.csv"
    db_file = "search_db/canonical_human_peptides_8_12.csv_mz_filtered.msp"
    charge = 1
    indistinguishability_from_db(
        search_pep_file,
        db_file,
        charge=charge,
        output_file=f"search_db/charlottes13peptides_canonical_human_peptides_8_12_charge_{charge}.csv",
        prefilter_mz=False,
    )
