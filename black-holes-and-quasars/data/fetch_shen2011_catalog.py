import pathlib

import requests

VIZIER_URL = "https://vizier.cds.unistra.fr/viz-bin/asu-tsv"

COLUMNS = [
    "SDSS", "RAJ2000", "DEJ2000", "z",
    "logLbol", "e_logLbol",
    "logL5100", "e_logL5100",
    "W(BHb)", "e_W(BHb)",
    "logBHHV", "logBHMM", "logBHCV",
    "logBH", "e_logBH",
    "BAL",
]

OUTPUT_PATH = pathlib.Path(__file__).parent / "shen2011_raw.tsv"


def fetch(row_limit: int = 30000) -> None:
    """Shen et al. 2011, ApJS 194, 45 (VizieR J/ApJS/194/45/catalog)."""
    params = {
        "-source": "J/ApJS/194/45/catalog",
        "-out.max": row_limit,
        "-out": ",".join(COLUMNS),
    }
    response = requests.get(VIZIER_URL, params=params, timeout=60)
    response.raise_for_status()
    OUTPUT_PATH.write_text(response.text)


if __name__ == "__main__":
    fetch()
