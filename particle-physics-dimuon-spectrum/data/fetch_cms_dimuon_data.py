import pathlib

import pandas as pd
import requests

SOURCE_URL = "https://raw.githubusercontent.com/cms-opendata-education/zboson-exercise/master/DoubleMuRun2011A.csv"
OUTPUT_PATH = pathlib.Path(__file__).parent / "dimuon_sample.csv"
SAMPLE_SIZE = 80000
RANDOM_SEED = 42


def fetch(sample_size: int = SAMPLE_SIZE) -> None:
    """CMS Open Data, Run2011A DoubleMu primary dataset (real 7 TeV LHC collisions)."""
    response = requests.get(SOURCE_URL, timeout=120)
    response.raise_for_status()
    full_path = pathlib.Path(__file__).parent / "_full_download.csv"
    full_path.write_bytes(response.content)

    df = pd.read_csv(full_path)
    sample = df.sample(n=sample_size, random_state=RANDOM_SEED).sort_values(["Run", "Event"])
    sample.to_csv(OUTPUT_PATH, index=False)
    full_path.unlink()


if __name__ == "__main__":
    fetch()
