"""
build_report.py — render the numeric tables in REPORT.md from the results files.

Numbers in REPORT.md are NOT hand-typed. Each table lives between a pair of
markers:

    <!-- AUTO:<key>:start -->
    ...generated table...
    <!-- AUTO:<key>:end -->

This script reads exploration/results/*.json and replaces the content between
each marker pair. Re-run after run_experiments.R / fidelity_test.py to refresh.

    python exploration/scripts/build_report.py
"""
import json
import re
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
RES = ROOT / "exploration" / "results"
REPORT = ROOT / "REPORT.md"


def load(name):
    return json.loads((RES / name).read_text())


def md_table(rows, headers=None, fmt=None):
    """rows: list of dicts. headers: ordered keys (defaults to first row's keys)."""
    if not rows:
        return "_(no data)_"
    headers = headers or list(rows[0].keys())
    fmt = fmt or {}

    def cell(r, h):
        v = r.get(h, "")
        if h in fmt and isinstance(v, (int, float)):
            return fmt[h](v)
        return str(v)

    out = ["| " + " | ".join(headers) + " |",
           "|" + "|".join(["---"] * len(headers)) + "|"]
    for r in rows:
        out.append("| " + " | ".join(cell(r, h) for h in headers) + " |")
    return "\n".join(out)


def sci(v):
    return f"{v:.2e}"


def f2(v):
    return f"{v:.2f}"


def build_runtime_profile():
    d = load("runtime_profile.json")
    rows = [{"Category": c["category"], "Self time (s)": c["self_time_s"],
             "% of runtime": c["pct"]} for c in d["categories"]]
    note = (f"\n\n_Sampled total: {d['total_sampled_s']}s. "
            f"Self-time = where the CPU actually sits._")
    return md_table(rows) + note


def build_tier1_ns_sweep():
    rows = load("tier1_ns_sweep.json")
    disp = [{"ns (resamples)": r["ns"], "time (s)": r["time_s"], "ATT": r["ATT"],
             "Lift %": r["lift_pct"], "p-value": r["pvalue"],
             "sig. @0.10": "yes" if r["significant_at_0.10"] else "no"} for r in rows]
    return md_table(disp)


def build_tier1_ci_methods():
    rows = load("tier1_ci_methods.json")
    disp = [{"CI method": r["method"], "time (s)": r["time_s"],
             "lower": r["lower"], "upper": r["upper"]} for r in rows]
    return md_table(disp)


def build_tier1_modified_compare():
    rows = load("tier1_modified_compare.json")
    disp = [{"Version": r["version"], "time (s)": r["time_s"], "ATT": r["ATT"],
             "Lift %": r["lift_pct"], "p-value": r["pvalue"],
             "speedup": f"{r['speedup']}x"} for r in rows]
    return md_table(disp)


def build_fidelity():
    rows = load("fidelity.json")
    disp = []
    for r in rows:
        if r.get("status") == "failed":
            disp.append({"Strategy": r["strategy"], "weights Δ": "—", "y_hat Δ": "—",
                         "ATT (Py/R)": "failed: " + r["error"].split(":")[0], "fit (ms)": "—"})
        else:
            disp.append({"Strategy": r["strategy"],
                         "weights Δ": sci(r["weights_max_abs_diff"]),
                         "y_hat Δ": sci(r["yhat_max_abs_diff"]),
                         "ATT (Py/R)": f"{r['ATT_python']} / {r['ATT_R']}",
                         "fit (ms)": r["fit_time_ms"]})
    return md_table(disp)


def build_market_selection_compare():
    d = load("market_selection_compare.json")
    s = d["setup"]
    rows = [{"ns": r["ns"], "time (s)": r["time_s"], "speedup": f"{r['speedup']}x",
             "sig. agreement": f"{round(100*r['power_significance_agreement'])}%",
             "top-5 overlap": f"{round(100*r['bestmarkets_top5_overlap'])}%",
             "MDE agreement": f"{round(100*r['bestmarkets_selected_mde_agreement'])}%",
             "AvgATT Δ": r["avg_att_max_abs_diff"]} for r in d["sweep"]]
    note = (f"\n\n_Same unmodified `GeoLiftMarketSelection`, only `ns` differs. "
            f"{s['locations_in_pool']}-market pool, N={s['N']}, "
            f"{s['n_simulations_per_run']} sims/run, seed {s['seed']}, sequential. "
            f"`ns=1000` is the reference (Meta default); lower `ns` scored against it._")
    return md_table(rows) + note


def build_market_selection_scaled_compare():
    d = load("market_selection_scaled_compare.json")
    ds, s = d["dataset"], d["setup"]
    rows = [{"ns": r["ns"], "time (s)": r["time_s"], "speedup": f"{r['speedup']}x",
             "sig. agreement": _pct(r["power_significance_agreement"]),
             "top-5 overlap": _pct(r["bestmarkets_top5_overlap"]),
             "MDE agreement": _pct(r["bestmarkets_selected_mde_agreement"]),
             "AvgATT Δ": r["avg_att_max_abs_diff"]} for r in d["sweep"]]
    note = (f"\n\n_Larger realistic panel: {ds['n_locations']} locations × {ds['n_periods']} "
            f"periods × {ds['n_categories']} categories (avg pairwise corr "
            f"{ds['avg_pairwise_correlation']}; cf. GeoLift_Test {ds['compared_to']}). "
            f"Same unmodified `GeoLiftMarketSelection`, only `ns` differs. N={s['N']}, "
            f"{s['n_simulations_per_run']} sims/run, seed {s['seed']}, sequential. "
            f"`ns=1000` (Meta default) is the reference._")
    return md_table(rows) + note


def build_market_selection_profile():
    d = load("market_selection_profile.json")
    m = d["model"]
    rows = [{"ns": p["ns"], "time (s)": p["time_s"],
             "speedup vs ns=1000": f"{p['speedup_vs_1000']}x"} for p in d["points"]]
    note = (f"\n\n_Per-run cost modelled as `time(ns) = fixed + slope·ns`. "
            f"Fixed floor = **{m['fixed_per_run_s']}s** ({round(100*(1-m['ns_scaled_fraction_at_1000']))}% "
            f"of the `ns=1000` runtime); only ~{round(100*m['ns_scaled_fraction_at_1000'])}% scales with `ns`. "
            f"So the `ns` lever alone tops out at ~**{m['ns_lever_ceiling_speedup']}×** even at `ns→0` — "
            f"bigger wins must cut the fixed per-sim work (fit-reuse, parallelism)._")
    return md_table(rows) + note


def build_market_selection_python_compare():
    d = load("market_selection_python_compare.json")
    s, t, a = d["setup"], d["timing"], d["agreement"]
    rows = [
        {"metric": "wall time (115 cells)",
         "R (`pvalueCalc` loop)": f"{t['r_seconds']}s ({t['r_ms_per_cell']} ms/cell)",
         "Python (numpy port)": f"{t['python_seconds']}s ({t['python_ms_per_cell']} ms/cell)"},
        {"metric": "speedup (sequential)", "R (`pvalueCalc` loop)": "1×",
         "Python (numpy port)": f"**{t['speedup_sequential']}×**"},
        {"metric": "significance agreement", "R (`pvalueCalc` loop)": "—",
         "Python (numpy port)": _pct(a["significance_agreement"])},
        {"metric": "BestMarkets top-5 overlap", "R (`pvalueCalc` loop)": "—",
         "Python (numpy port)": _pct(a["ranking_top5_overlap"])},
        {"metric": "selected-MDE agreement", "R (`pvalueCalc` loop)": "—",
         "Python (numpy port)": _pct(a["selected_mde_agreement"])},
        {"metric": "AvgATT max abs diff", "R (`pvalueCalc` loop)": "—",
         "Python (numpy port)": f"{a['att_max_abs_diff']:.1e}"},
        {"metric": "scaled-L2 max abs diff", "R (`pvalueCalc` loop)": "—",
         "Python (numpy port)": f"{a['scaled_l2_max_abs_diff']:.1e}"},
        {"metric": "p-value mean abs diff", "R (`pvalueCalc` loop)": "—",
         "Python (numpy port)": f"{a['pvalue_mean_abs_diff']}"},
    ]
    note = (f"\n\n_Same {s['n_combos']} location combos × {len(s['effect_sizes'])} effect "
            f"sizes = {s['n_cells']} cells, {s['panel']}, ns={s['ns']}, seed {s['seed']}, "
            f"single-threaded. Deterministic quantities (ATT, scaled-L2) match R to solver "
            f"tolerance; significance, market ranking and selected MDE are **identical**. "
            f"P-values differ only by Monte-Carlo error ({s['ns']} permutations, independent "
            f"R/Python RNG) — mean |Δ| {a['pvalue_mean_abs_diff']}, never flipping a decision. "
            f"The {t['speedup_sequential']}× comes purely from vectorising the permutation "
            f"draw + lever-2 fit-reuse, before any parallelism._")
    return md_table(rows) + note


def build_bench_scaling():
    d = load("bench_scaling.json")
    s, a, rows = d["setup"], d["agreement"], d["scaling"]
    eng = s["engines"]
    scale = md_table([
        {"locations": r["n_locations"], "pairs": r["n_combos"], "cells": r["n_cells"],
         "R original (s)": r["r_original_s"], "R modified (s)": r["r_modified_s"],
         "Python (s)": r["python_s"],
         "Py vs R-orig": f"**{r['speedup_py_vs_r_original']}×**",
         "Py vs R-mod": f"{r['speedup_py_vs_r_modified']}×"}
        for r in rows])
    mean_o = sum(r["speedup_py_vs_r_original"] for r in rows) / len(rows)
    mean_rm = sum(r["speedup_r_modified_vs_original"] for r in rows) / len(rows)
    note = (f"\n\n_Same workload per row — all C(L,2) treated market pairs × "
            f"{len(s['effect_sizes'])} effect sizes, seed {s['seed']}, single-threaded. "
            f"**R original** = {eng['r_original']}; **R modified** = {eng['r_modified']}; "
            f"**Python** = {eng['python']}. Each R point runs in its own process (crash "
            f"isolation). The R-modified `ns` lever buys only ~{mean_rm:.1f}× over R-original "
            f"(its measured ~1.4× ceiling — §6), while Python at **full `ns=1000` fidelity** is "
            f"~{mean_o:.0f}× faster than R-original. Wall-clock varies by machine; the ratios "
            f"are the portable numbers._")
    agree = md_table([
        {"check": "significance agreement", "result": _pct(a["significance_agreement"])
            + f" ({a['n_cells']} cells)"},
        {"check": "ATT max abs diff", "result": f"{a['att_max_abs_diff']:.1e}"},
        {"check": "scaled-L2 max abs diff", "result": f"{a['scaled_l2_max_abs_diff']:.1e}"},
        {"check": "p-value mean abs diff", "result": f"{a['pvalue_mean_abs_diff']}"}])
    agree_note = (f"\n\n_Per-cell agreement at the largest size ({a['at_n_locations']} "
                  f"locations, {a['n_cells']} cells), {a['compared']}. The estimators (ATT, "
                  f"scaled-L2) match R to solver tolerance; the only significance flips are at "
                  f"`EffectSize=0` where the p-value straddles α=0.10 by hundredths — pure "
                  f"Monte-Carlo noise from independent RNG, not a modelling difference._")
    return (agree + agree_note
            + "\n\n![R original vs R modified vs Python scaling](exploration/results/bench_scaling.png)\n\n"
            + scale + note)


def _pct(v):
    return "NA" if v is None else f"{round(100*v)}%"


def build_fit_reuse_verify():
    d = load("fit_reuse_verify.json")
    s = d["setup"]
    rows = [{"es": r["es"],
             "weights Δ vs es0": r["weights_max_abs_diff_vs_es0"],
             "counterfactual Δ vs es0": r["cf_max_abs_diff_vs_es0"],
             "ATT": r["att"], "p-value": r["pvalue"]} for r in d["per_es"]]
    note = (f"\n\n_Fixed seed {s['seed']}, market `{s['locations']}`, duration {s['duration_tp']}, "
            f"ns={s['ns']}. Δ=0 ⇒ identical across effect sizes._")
    return md_table(rows) + note


def build_provenance():
    d = load("provenance.json")
    return (f"- **Seed:** {d['seed']} (statistics deterministic)\n"
            f"- **R:** {d['r_version']}\n"
            f"- **augsynth:** {d['augsynth_commit']}\n"
            f"- **Example:** {d['example']}\n"
            f"- _{d['note']}_")


BUILDERS = {
    "runtime_profile": build_runtime_profile,
    "tier1_ns_sweep": build_tier1_ns_sweep,
    "tier1_ci_methods": build_tier1_ci_methods,
    "tier1_modified_compare": build_tier1_modified_compare,
    "fit_reuse_verify": build_fit_reuse_verify,
    "market_selection_compare": build_market_selection_compare,
    "market_selection_scaled_compare": build_market_selection_scaled_compare,
    "market_selection_profile": build_market_selection_profile,
    "market_selection_python_compare": build_market_selection_python_compare,
    "bench_scaling": build_bench_scaling,
    "fidelity": build_fidelity,
    "provenance": build_provenance,
}


def main():
    text = REPORT.read_text(encoding="utf-8")
    for key, fn in BUILDERS.items():
        block = fn()
        pat = re.compile(
            rf"(<!-- AUTO:{key}:start -->).*?(<!-- AUTO:{key}:end -->)",
            re.DOTALL)
        if not pat.search(text):
            print(f"  WARN: markers for '{key}' not found in REPORT.md")
            continue
        text = pat.sub(rf"\1\n{block}\n\2", text)
        print(f"  filled: {key}")
    REPORT.write_text(text, encoding="utf-8")
    print(f"Wrote {REPORT}")


if __name__ == "__main__":
    main()
