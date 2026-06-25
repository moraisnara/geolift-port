"""
_compare_common.py — shared paths + tiny helpers for the geolift_py comparison /
benchmark scripts (compare_market_selection.py, compare_city_example.py,
validate_market_selection.py, bench_scaling.py).

This is NOT part of the installable `geolift_fast` package — it only factors out
boilerplate the analysis scripts repeated verbatim. Importable when a script is
run as `python geolift_py/<script>.py` (the script's dir is on sys.path).
"""
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
RES = ROOT / "exploration" / "results"
DATA = ROOT / "exploration" / "data"
BENCH = RES / "bench"


def sig3(x):
    """3 significant figures as a JSON-friendly float — the `float(f"{x:.3g}")`
    idiom the compare scripts use for reporting max-abs-diffs."""
    return float(f"{x:.3g}")


def read_lower(path):
    """Read a CSV and lowercase its `location` column, matching Panel /
    GeoDataRead lowercasing so R tables join cleanly to Python output."""
    df = pd.read_csv(path)
    df["location"] = df["location"].str.lower()
    return df
