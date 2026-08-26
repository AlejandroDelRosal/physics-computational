import pathlib

import pandas as pd

RAW_PATH = pathlib.Path(__file__).parent.parent / "data" / "shen2011_raw.tsv"

NUMERIC_COLUMNS = [
    "RAJ2000", "DEJ2000", "z",
    "logLbol", "e_logLbol",
    "logL5100", "e_logL5100",
    "W(BHb)", "e_W(BHb)",
    "logBHHV", "logBHMM", "logBHCV",
    "logBH", "e_logBH",
    "BAL",
]


def load_raw_catalog(path: pathlib.Path = RAW_PATH) -> pd.DataFrame:
    df = pd.read_csv(path, sep="\t", comment="#", skiprows=[1, 2])
    df.columns = [c.strip() for c in df.columns]
    for col in NUMERIC_COLUMNS:
        df[col] = pd.to_numeric(df[col].astype(str).str.strip(), errors="coerce")
    df["SDSS"] = df["SDSS"].astype(str).str.strip()
    return df


def load_hbeta_subsample(path: pathlib.Path = RAW_PATH) -> pd.DataFrame:
    df = load_raw_catalog(path)
    has_hbeta = df["W(BHb)"].notna() & df["logL5100"].notna() & df["logBHHV"].notna()
    not_bal = df["BAL"] == 0
    return df.loc[has_hbeta & not_bal].reset_index(drop=True)
