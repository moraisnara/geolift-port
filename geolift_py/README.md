# geolift-fast

A fast Python port of [GeoLift](https://github.com/facebookincubator/GeoLift), covering
**both halves** of the workflow: **market-selection / power inference** (pre-test) and the
**single-test lift measurement** `GeoLift()` (post-test).

It reproduces GeoLift's per-cell decisions — ATT, scaled-L2 imbalance, the conformal
"iid" permutation p-value, and the selected-market ranking — to solver tolerance
(~1e-8 on the estimators, significance/ranking identical), while running **~10–14×
faster than the R inner loop single-threaded**. The speedup is structural, not a
parameter trade-off: the `ns` conformal permutations are drawn as one numpy matrix
instead of R's interpreted `sapply(1:ns, …)`, and the effect-size-invariant SCM fit is
computed once per market combo and reused across every effect size ("lever 2").

The post-test `geolift()` mirrors R's `GeoLift()` on an augsynth-faithful engine
(fixed-effects + ridge SCM, conformal p-value, conformal-CI grid inversion with jackknife+
fallback): deterministic quantities match R to solver tolerance, the conformal p-value
matches in the limit (identical residual vector), ~12.6× faster on the CI task.

See [`../REPORT.md`](../REPORT.md) §7 (pre-test) and §7d (post-test) for validation and benchmarks.

## Install

From GitHub (e.g. on a GCP VM), pinned to a commit/tag for reproducibility:

```bash
pip install "git+https://github.com/moraisnara/geolift-port.git@v0.1.0"
```

The package lives in a subdirectory of this repo, so `pip` is pointed at it via the
`pyproject.toml` at the repo root — the command above works as-is. Dependencies
(`numpy`, `pandas`, `scipy`, `osqp`) install automatically; add the plotting extra with
`pip install "git+…#egg=geolift-fast[plot]"` if you want matplotlib.

Local editable install for development:

```bash
git clone https://github.com/moraisnara/geolift-port.git && cd geolift-port
pip install -e .
```

### Minimal GCP recipe

```bash
sudo apt-get update && sudo apt-get install -y python3-venv git
python3 -m venv ~/venv && source ~/venv/bin/activate
pip install "git+https://github.com/moraisnara/geolift-port.git@v0.1.0"
python my_market_selection.py
```

## Usage

```python
import pandas as pd
from geolift_fast import power_curves, best_markets, all_pairs, Panel

# Long/tidy panel: columns location, date, Y (one row per location-period).
df = pd.read_csv("panel.csv")
panel = Panel.from_long_df(df)

# Candidate treated markets: all pairs, or pass your own list of lists.
combos = all_pairs(panel, size=2)                 # e.g. [["chicago","portland"], ...]

# Run the power simulation -> tidy PowerCurves DataFrame.
pc = power_curves(
    panel, combos,
    treatment_periods=14,                         # int or list of test-window lengths
    effect_sizes=[-0.1, -0.05, 0.0, 0.05, 0.1],
    ns=1000,                                       # conformal permutations (Meta default)
    alpha=0.10,
    seed=42,                                       # statistics are deterministic given seed
)
# columns: location, duration, EffectSize, pvalue, power, AvgATT,
#          AvgDetectedLift, AvgScaledL2Imbalance

# Rank markets by their MDE: the significant effect closest to zero in either
# direction (can be negative), then a composite rank a la R's BestMarkets.
ranking = best_markets(pc, alpha=0.10)            # columns: rank, location,
#          duration, MDE, AvgDetectedLift, abs_lift_in_zero
print(ranking.head())
```

### Lower-level building blocks

```python
from geolift_fast import Panel, ComboFit, simulate_combo, scm_weights
import numpy as np

panel = Panel.from_long_csv("panel.csv")
rng = np.random.default_rng(42)
cells = simulate_combo(panel, ["chicago", "portland"], tp=14,
                       effect_sizes=[0.0, 0.05, 0.1], ns=1000, rng=rng, alpha=0.1)
```

### Post-test measurement (`GeoLift()` equivalent)

```python
from geolift_fast import Panel, geolift

panel = Panel.from_long_csv("panel.csv")
r = geolift(
    panel, locations=["chicago", "portland"],
    treatment_start_time=91, treatment_end_time=105,   # 1-based period indices
    model="none",                                       # or "ridge" / "best"
    confidence_intervals=True, method="conformal",      # falls back to jackknife+ if conformal fails
    ns=1000, alpha=0.10, seed=42,
)
print(r.summary())          # ATT, Percent Lift, Incremental, p-value, scaled-L2, CI
r.att, r.percent_lift, r.pvalue, r.lower_conf_int, r.upper_conf_int, r.y_hat
```

## Fidelity & scope

- **Pre-test config covered:** the GeoLift market-selection default — `fixed_effects=TRUE`,
  `model="none"`, `conformal_type="iid"`, two-sided test (`stat = sum|x|`).
- **Post-test config covered:** `model` in {`none`, `ridge`, `best`}, `fixed_effects=TRUE`,
  conformal "iid" p-value, conformal-CI grid inversion + jackknife+ fallback, stat tests
  Total/Negative/Positive.
- **Validated against R** to ~1e-11 on deterministic pre-test internals and to solver
  tolerance on the post-test estimators (ATT, weights, counterfactual, scaled-L2, ridge λ,
  jackknife+ CI). Significance and ranking identical; conformal p-values differ only by
  Monte-Carlo error (independent RNG, identical limiting value), never flipping a decision
  except at the α boundary under the null.
- **Out of scope:** `model="GSYN"` (augsynth delegates it to the separate `gsynth` factor-model
  package), non-iid conformal variants, and plotting helpers.
