# geolift-fast

A fast Python port of [GeoLift](https://github.com/facebookincubator/GeoLift)'s
**market-selection / power inference** inner loop.

It reproduces GeoLift's per-cell decisions — ATT, scaled-L2 imbalance, the conformal
"iid" permutation p-value, and the selected-market ranking — to solver tolerance
(~1e-8 on the estimators, significance/ranking identical), while running **~10–14×
faster than the R inner loop single-threaded**. The speedup is structural, not a
parameter trade-off: the `ns` conformal permutations are drawn as one numpy matrix
instead of R's interpreted `sapply(1:ns, …)`, and the effect-size-invariant SCM fit is
computed once per market combo and reused across every effect size ("lever 2").

See [`../REPORT.md`](../REPORT.md) §7 for the validation and the R-vs-Python scaling benchmark.

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
# columns: location, duration, EffectSize, pvalue, power, AvgATT, AvgScaledL2Imbalance

# Rank markets by minimum detectable effect (smallest significant positive lift).
ranking = best_markets(pc, alpha=0.10)            # columns: rank, location, duration, MDE
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

## Fidelity & scope

- **Config covered:** the GeoLift market-selection default — `fixed_effects=TRUE`,
  `model="none"`, `conformal_type="iid"`, two-sided test (`stat = sum|x|`).
- **Validated against R** to ~1e-11 on deterministic internals (weights, residuals,
  observed statistic) and ~1e-8 on ATT; significance and ranking identical. P-values
  differ only by Monte-Carlo error (independent RNG), never flipping a decision except
  at the α boundary under the null.
- The `augsynth` augmentation (`progfunc != "none"`) and non-iid conformal variants are
  out of scope for this port.
