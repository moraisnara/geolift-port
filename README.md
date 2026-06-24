# geolift-port

Porting and optimizing Meta's [GeoLift](https://github.com/facebookincubator/GeoLift)
market-selection / power inference from R to Python — and measuring, fixest-style,
exactly how much faster the port is and whether it makes the **same decisions**.

The headline: on the market-selection inner loop, the Python port runs **~10–11×
faster than the unmodified R engine single-threaded**, while reproducing GeoLift's
per-cell ATT, scaled-L2 imbalance, conformal p-value, and selected-market ranking to
solver tolerance. The speedup is **structural** (vectorized conformal permutations +
effect-size-invariant fit reuse), not a quality trade-off — see
[`REPORT.md`](REPORT.md) §7 for the validation and the scaling benchmark.

![scaling benchmark](exploration/results/bench_scaling.png)

## Three engines compared

| Folder | Engine | What it is |
|---|---|---|
| [`geolift_r_original/`](geolift_r_original/) | **R original** | Meta's GeoLift v2.7.5, unmodified (trimmed to runtime essentials), `ns=1000` Meta default |
| [`geolift_r_modified/`](geolift_r_modified/) | **R modified** | Parameter presets (`fast_presets.R`) passed to the *unchanged* engine — e.g. `ns=100`. Not a fork |
| [`geolift_py/`](geolift_py/) | **Python port** | `geolift_fast`, an installable Python library reproducing the inner loop, `ns=1000` full fidelity |

The R-modified preset only buys ~1.5× (the `ns` lever is bounded — ~70% of per-sim
cost is fixed), while the Python port sits a full structural tier below both. That
contrast is the point of the benchmark.

## Repo layout

```
geolift_r_original/   trimmed vendored copy of Meta GeoLift v2.7.5 (baseline)
geolift_r_modified/   fast_presets.R — parameter presets, not a code fork
geolift_py/           the Python port + the installable `geolift_fast` library
exploration/          reproducibility harness:
  data/                 fixed input panels (seeded)
  scripts/              R + Python scripts that regenerate every result
  results/              generated JSON/CSV/PNG the report injects from
  ground_truth/         frozen R ground-truth cells
  logs/                 run transcripts (provenance)
pyproject.toml        makes geolift_py/geolift_fast pip-installable
REPORT.md             full write-up: method, fidelity, benchmarks, reproduce steps
```

## Quickstart

**Python library** (works on a GCP VM):

```bash
pip install "git+https://github.com/moraisnara/geolift-port.git@v0.1.0"
```

```python
from geolift_fast import power_curves, best_markets, all_pairs, Panel

panel = Panel.from_long_csv("panel.csv")          # columns: location, date, Y
pc = power_curves(panel, all_pairs(panel, size=2),
                  treatment_periods=14, effect_sizes=[-0.1, 0, 0.1], ns=1000, seed=42)
ranking = best_markets(pc, alpha=0.10)            # rank, location, duration, MDE
```

**Reproduce the benchmark** (numbers in `REPORT.md` are injected from
`exploration/results/`, never hand-typed; everything is seeded with `SEED=42`):

```bash
# R points (one fresh process per size; ns=1000 = original, ns=100 = modified):
Rscript exploration/scripts/bench_scaling_R.R 14 1000
Rscript exploration/scripts/bench_scaling_R.R 14 100
# Python timings + agreement + plot:
python geolift_py/bench_scaling.py
# Re-render REPORT.md from the results:
python exploration/scripts/build_report.py
```

See [`geolift_py/README.md`](geolift_py/README.md) for the full library API and a GCP
install recipe, and [`REPORT.md`](REPORT.md) for the complete study.

## License

MIT — see [`LICENSE`](LICENSE). The vendored Meta GeoLift code under
`geolift_r_original/` retains its own MIT license (Copyright Meta Platforms, Inc.).
