import json
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy import stats


ROOT = Path(__file__).resolve().parents[1] / "results" / "raw_runs"
FIG = ROOT / "figures"
FIG.mkdir(exist_ok=True)
RNG = np.random.default_rng(20261008)

ALGORITHMS = ["CPO-DBO", "CPO", "DBO", "GA", "IPSO", "GA-PSO", "ACO-GA"]
ABLATIONS = [
    "No 2-opt", "Fixed probabilities", "No group bias", "No elitism",
    "Swap-only mutation", "No personal-best guidance",
]


def holm(pvalues):
    pvalues = np.asarray(pvalues, float)
    order = np.argsort(pvalues)
    adjusted = np.empty_like(pvalues)
    running = 0.0
    m = len(pvalues)
    for rank, idx in enumerate(order):
        running = max(running, (m - rank) * pvalues[idx])
        adjusted[idx] = min(1.0, running)
    return adjusted


def rank_biserial(diff):
    diff = np.asarray(diff)
    nonzero = diff[np.abs(diff) > 1e-14]
    ranks = stats.rankdata(np.abs(nonzero))
    denom = ranks.sum()
    return float((ranks[nonzero > 0].sum() - ranks[nonzero < 0].sum()) / denom) if denom else 0.0


def bootstrap_mean_ci(values, n_boot=50_000):
    values = np.asarray(values, float)
    chunks = []
    remaining = n_boot
    while remaining:
        take = min(5000, remaining)
        sample = RNG.choice(values, size=(take, values.size), replace=True).mean(axis=1)
        chunks.append(sample)
        remaining -= take
    boot = np.concatenate(chunks)
    return [float(x) for x in np.quantile(boot, [0.025, 0.975])]


def exact_sign_flip_p(diff):
    diff = np.asarray(diff, float)
    observed = abs(diff.mean())
    n = diff.size
    extreme = 0
    total = 1 << n
    for start in range(0, total, 65536):
        ids = np.arange(start, min(start + 65536, total), dtype=np.uint32)
        sums = np.zeros(ids.size)
        for bit in range(n):
            signs = 1.0 - 2.0 * ((ids >> bit) & 1)
            sums += signs * diff[bit]
        extreme += int(np.sum(np.abs(sums / n) >= observed - 1e-14))
    return extreme / total


def pair_result(instance_means, reference, comparator):
    pivot = instance_means.pivot(index="instance", columns="variant", values="objective_au")
    diff = pivot[comparator].to_numpy() - pivot[reference].to_numpy()
    wil = stats.wilcoxon(diff, alternative="two-sided", method="auto")
    return {
        "reference": reference,
        "comparator": comparator,
        "mean_difference_comparator_minus_reference": float(diff.mean()),
        "bootstrap_95_ci": bootstrap_mean_ci(diff),
        "wilcoxon_statistic": float(wil.statistic),
        "wilcoxon_p": float(wil.pvalue),
        "rank_biserial": rank_biserial(diff),
        "reference_wins": int(np.sum(diff > 0)),
        "ties": int(np.sum(np.isclose(diff, 0))),
        "reference_losses": int(np.sum(diff < 0)),
    }


def read_curves(variants):
    rows = []
    for variant in variants:
        label = (variant.replace(" ", "_").replace("=", "-").replace("%", "pct")
                 .replace("/", "_").replace("–", "-").replace("-", "_"))
        for path in sorted((ROOT / "runs" / label).glob("*.npz")):
            with np.load(path) as data:
                rows.append((variant, data["convergence_progress"].astype(float)))
    return rows


def save_figure(fig, name):
    fig.savefig(FIG / f"{name}.png", dpi=400, bbox_inches="tight")
    fig.savefig(FIG / f"{name}.svg", bbox_inches="tight")
    fig.savefig(FIG / f"{name}.pdf", bbox_inches="tight")
    plt.close(fig)


def main():
    runs = pd.read_csv(ROOT / "all_study_runs.csv")
    assert len(runs) == 4000
    assert np.allclose(runs["objective_au"], runs["recomputed_objective_au"], atol=1e-9)
    assert runs.groupby(["variant", "instance"]).size().eq(10).all()
    instance_means = runs.groupby(["study", "variant", "instance"], as_index=False).agg(
        objective_au=("objective_au", "mean"),
        objective_sd=("objective_au", "std"),
        seconds=("seconds", "mean"),
        first_hit_iteration=("first_hit_iteration", "median"),
    )
    instance_means.to_csv(ROOT / "instance_level_means.csv", index=False, encoding="utf-8-sig")

    summary = runs.groupby(["study", "variant"], as_index=False).agg(
        mean_objective_au=("objective_au", "mean"),
        sd_objective_au=("objective_au", "std"),
        median_objective_au=("objective_au", "median"),
        mean_seconds=("seconds", "mean"),
        median_first_hit=("first_hit_iteration", "median"),
        evaluations=("full_tour_evaluations", "median"),
    )
    summary.to_csv(ROOT / "study_summary.csv", index=False, encoding="utf-8-sig")

    comp = instance_means[instance_means.variant.isin(ALGORITHMS)]
    comp_pivot = comp.pivot(index="instance", columns="variant", values="objective_au")[ALGORITHMS]
    friedman = stats.friedmanchisquare(*[comp_pivot[c].to_numpy() for c in ALGORITHMS])
    comp_pairs = [pair_result(comp, "CPO-DBO", name) for name in ALGORITHMS[1:]]
    adjusted = holm([row["wilcoxon_p"] for row in comp_pairs])
    for row, value in zip(comp_pairs, adjusted):
        row["holm_p"] = float(value)
    pd.DataFrame(comp_pairs).to_csv(ROOT / "algorithm_pairwise_tests.csv", index=False, encoding="utf-8-sig")

    abl_variants = ["CPO-DBO"] + ABLATIONS
    abl = instance_means[instance_means.variant.isin(abl_variants)]
    abl_pairs = [pair_result(abl, "CPO-DBO", name) for name in ABLATIONS]
    for row in abl_pairs:
        pivot = abl.pivot(index="instance", columns="variant", values="objective_au")
        diff = pivot[row["comparator"]].to_numpy() - pivot["CPO-DBO"].to_numpy()
        row["exact_sign_flip_p"] = float(exact_sign_flip_p(diff))
    adjusted = holm([row["exact_sign_flip_p"] for row in abl_pairs])
    for row, value in zip(abl_pairs, adjusted):
        row["holm_sign_flip_p"] = float(value)
    pd.DataFrame(abl_pairs).to_csv(ROOT / "ablation_pairwise_tests.csv", index=False, encoding="utf-8-sig")

    sensitivity = {
        "Population size": ["M=200", "CPO-DBO", "M=600"],
        "Schedule scale": ["Tscale=5000", "Tscale=15000", "CPO-DBO"],
        "Elite ratio": ["Elite=1%", "CPO-DBO", "Elite=5%"],
        "2-opt attempts": ["No 2-opt", "CPO-DBO", "2-opt=3"],
    }
    sensitivity_rows = []
    sensitivity_tests = []
    sensitivity_pairs = []
    for factor, levels in sensitivity.items():
        subset = instance_means[instance_means.variant.isin(levels)]
        pivot = subset.pivot(index="instance", columns="variant", values="objective_au")[levels]
        test = stats.friedmanchisquare(*[pivot[c].to_numpy() for c in levels])
        sensitivity_tests.append({"factor": factor, "friedman_chi2": float(test.statistic), "p": float(test.pvalue)})
        factor_pairs = [pair_result(subset, "CPO-DBO", level) for level in levels if level != "CPO-DBO"]
        factor_adjusted = holm([row["wilcoxon_p"] for row in factor_pairs])
        for row, value in zip(factor_pairs, factor_adjusted):
            row["factor"] = factor
            row["holm_p_within_factor"] = float(value)
            sensitivity_pairs.append(row)
        for level in levels:
            values = pivot[level].to_numpy()
            sensitivity_rows.append({
                "factor": factor,
                "level": level,
                "mean_objective_au": float(values.mean()),
                "sd_across_instances": float(values.std(ddof=1)),
                "bootstrap_95_ci": json.dumps(bootstrap_mean_ci(values)),
            })
    pd.DataFrame(sensitivity_rows).to_csv(ROOT / "sensitivity_summary.csv", index=False, encoding="utf-8-sig")
    pd.DataFrame(sensitivity_tests).to_csv(ROOT / "sensitivity_friedman_tests.csv", index=False, encoding="utf-8-sig")
    pd.DataFrame(sensitivity_pairs).to_csv(ROOT / "sensitivity_pairwise_tests.csv", index=False, encoding="utf-8-sig")

    ranks = comp_pivot.rank(axis=1, method="average", ascending=True).mean(axis=0).sort_values()
    statistics = {
        "comparison": {
            "friedman_chi2": float(friedman.statistic),
            "friedman_p": float(friedman.pvalue),
            "average_ranks": {k: float(v) for k, v in ranks.items()},
            "pairwise_vs_cpo_dbo": comp_pairs,
        },
        "ablation": abl_pairs,
        "sensitivity": sensitivity_tests,
        "sensitivity_pairwise": sensitivity_pairs,
    }
    (ROOT / "statistics.json").write_text(json.dumps(statistics, ensure_ascii=False, indent=2), encoding="utf-8")

    # ---- Global style: Times New Roman, no grid, no minor ticks, ticks inward ----
    plt.rcParams.update({
        "font.family": "serif",
        "font.serif": ["Times New Roman", "Nimbus Roman", "DejaVu Serif"],
        "mathtext.fontset": "stix",
        "font.size": 9,
        "axes.linewidth": 0.8,
        "axes.grid": False,
       # "axes.xtick": False,
        "xtick.direction": "in",
        "ytick.direction": "in",
        "xtick.minor.visible": False,
        "ytick.minor.visible": False,
        "xtick.major.size": 3.5,
        "ytick.major.size": 3.5,
        "xtick.top": False,
        "ytick.right": False,
    })

    colors = ["#0072B2" if x == "CPO-DBO" else "#B8BEC7" for x in ALGORITHMS]
    fig, ax = plt.subplots(figsize=(7.2, 3.8))
    values = [comp_pivot[x].to_numpy() for x in ALGORITHMS]
    box = ax.boxplot(values, labels=ALGORITHMS, patch_artist=True, showfliers=False)
    for patch, color in zip(box["boxes"], colors):
        patch.set_facecolor(color)
        patch.set_alpha(0.85)
    for i, vals in enumerate(values, 1):
        jitter = RNG.normal(0, 0.045, size=len(vals))
        ax.scatter(np.full(len(vals), i) + jitter, vals, s=12, color="#222222", alpha=0.55, zorder=3)
    ax.set_ylabel("Relative exposure cost (a.u.)")
    ax.set_xlabel("Algorithm")
    # grid removed
    fig.tight_layout()
    save_figure(fig, "Fig19_algorithm_comparison")

    curves = read_curves(ALGORITHMS)
    fig, ax = plt.subplots(figsize=(7.2, 3.8))
    palette = dict(zip(ALGORITHMS, ["#0072B2", "#D55E00", "#009E73", "#CC79A7", "#E69F00", "#56B4E9", "#6A3D9A"]))
    x = np.arange(101)
    for algorithm in ALGORITHMS:
        matrix = np.vstack([curve for name, curve in curves if name == algorithm])
        ax.plot(x, matrix.mean(axis=0), label=algorithm, lw=1.6, color=palette[algorithm])
    ax.set_xlabel("Objective-evaluation budget used (%)")
    ax.set_ylabel("Best relative exposure cost (a.u.)")
    # grid removed
    ax.legend(ncol=4, frameon=False, fontsize=8, loc="upper center",
              bbox_to_anchor=(0.5, -0.13))
    ax.set_ylim(335, 445)
    inset = ax.inset_axes([0.53, 0.38, 0.44, 0.55])
    for algorithm in ALGORITHMS:
        matrix = np.vstack([curve for name, curve in curves if name == algorithm])
        inset.plot(x, matrix.mean(axis=0), lw=1.0, color=palette[algorithm])
    inset.set_xlim(0, 12)
    inset.set_ylim(300, 1850)
    inset.set_title("Initial phase", fontsize=8)
    inset.tick_params(labelsize=7)
    # inset grid removed
    fig.tight_layout()
    save_figure(fig, "Fig20_budget_normalized_convergence")

    abl_df = pd.DataFrame(abl_pairs).sort_values("mean_difference_comparator_minus_reference")
    fig, ax = plt.subplots(figsize=(7.2, 3.8))
    y = np.arange(len(abl_df))
    means = abl_df["mean_difference_comparator_minus_reference"].to_numpy()
    cis = np.vstack(abl_df["bootstrap_95_ci"].to_numpy())
    ax.errorbar(means, y, xerr=np.vstack([means - cis[:, 0], cis[:, 1] - means]), fmt="o", color="#0072B2", capsize=3)
    ax.axvline(0, color="#444444", lw=1, ls="--")
    ax.set_yticks(y, abl_df["comparator"])
    ax.set_xlabel("Cost increase after ablation (a.u.)")
    # grid removed
    fig.tight_layout()
    save_figure(fig, "Fig21_ablation_effects")

    sens_df = pd.DataFrame(sensitivity_rows)
    fig, axes = plt.subplots(2, 2, figsize=(7.2, 5.8))
    for ax, (factor, levels) in zip(axes.flat, sensitivity.items()):
        part = sens_df[sens_df.factor == factor].set_index("level").loc[levels]
        means = part.mean_objective_au.to_numpy()
        cis = np.vstack([json.loads(x) for x in part.bootstrap_95_ci])
        ax.errorbar(np.arange(3), means, yerr=np.vstack([means - cis[:, 0], cis[:, 1] - means]),
                    marker="o", color="#0072B2", capsize=3)
        labels = levels.copy()
        labels = ["M=400" if x == "CPO-DBO" and factor == "Population size" else x for x in labels]
        labels = ["Tscale=30000" if x == "CPO-DBO" and factor == "Schedule scale" else x for x in labels]
        labels = ["Elite=2%" if x == "CPO-DBO" and factor == "Elite ratio" else x for x in labels]
        labels = ["2-opt=1" if x == "CPO-DBO" and factor == "2-opt attempts" else x for x in labels]
        ax.set_xticks(np.arange(3), labels, rotation=12)
        ax.set_title(factor)
        ax.set_ylabel("Cost (a.u.)")
        # grid removed
    fig.tight_layout()
    save_figure(fig, "Fig22_parameter_sensitivity")

    checks = {
        "run_count": int(len(runs)),
        "expected_run_count": 4000,
        "objective_recalculation_max_error": float(np.max(np.abs(runs.objective_au - runs.recomputed_objective_au))),
        "all_variant_instance_counts_are_10": bool(runs.groupby(["variant", "instance"]).size().eq(10).all()),
        "figures": sorted(p.name for p in FIG.glob("*.png")),
    }
    (ROOT / "analysis_checks.json").write_text(json.dumps(checks, ensure_ascii=False, indent=2), encoding="utf-8")
    print(json.dumps(statistics, ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()