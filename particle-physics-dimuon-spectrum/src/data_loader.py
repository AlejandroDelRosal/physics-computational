import pathlib

import pandas as pd

DATA_PATH = pathlib.Path(__file__).parent.parent / "data" / "dimuon_sample.csv"


def load_dimuon_events(path: pathlib.Path = DATA_PATH) -> pd.DataFrame:
    df = pd.read_csv(path)
    return df[df["Q1"] * df["Q2"] < 0].reset_index(drop=True)
