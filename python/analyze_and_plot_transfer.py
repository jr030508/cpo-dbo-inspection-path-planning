import csv
import json
from pathlib import Path

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
from PIL import Image
from scipy.stats import friedmanchisquare, rankdata, wilcoxon

RAW = Path(__file__).resolve().parents[1] / "results" / "raw_runs"
ROOT = Path(__file__).resolve().parents[1] / "output" / "analyze_and_plot_transfer"
ROOT.mkdir(parents=True, exist_ok=True)
OLD = Path(__file__).resolve().parents[1] / "results" / "raw_runs"
FIG = ROOT / "figures"
FIG.mkdir(exist_ok=True)

BLUE = "#0072B2"
ORANGE = "#E69F00"
GREEN = "#009E73"
RED = "#D55E00"
GRAY = "#666666"
LIGHT_BLUE = "#56B4E9"


def apply_style():
    mpl.rcParams.update({
        "font.family": "serif",
        "font.serif": ["Times New Roman"],
        "mathtext.fontset": "stix",
        "font.size": 9,
        "axes.titlesize": 10.5,
        "axes.labelsize": 10,
        "xtick.labelsize": 8.5,
        "ytick.labelsize": 8.5,
        "legend.fontsize": 8.5,
        "axes.linewidth": 0.9,
        "figure.facecolor": "white",
        "axes.facecolor": "white",
        "savefig.facecolor": "white",
        "savefig.transparent": False,
        "axes.grid": False,
        "pdf.fonttype": 42,
        "ps.fonttype": 42,
    })


def style_axis(ax):
    ax.grid(False)
    for side in ("top", "right", "bottom", "left"):
        ax.spines[side].set_visible(True)
        ax.spines[side].set_linewidth(0.9)
        ax.spines[side].set_color("black")
    ax.minorticks_off()
    ax.tick_params(axis="both", which="minor", bottom=False, left=False, top=False, right=False)
    ax.tick_params(axis="both", which="major", direction="in", length=4.0, width=0.9,
                   bottom=True, left=True, top=False, right=False, pad=4)
    # Keep the top and right frame but remove their tick marks and labels.
    ax.tick_params(axis="x", top=False, labeltop=False)
    ax.tick_params(axis="y", right=False, labelright=False)


def save_all(fig, name):
    fig.savefig(FIG / f"{name}.png", dpi=600, bbox_inches="tight", pad_inches=0.03)
    fig.savefig(FIG / f"{name}.tiff", dpi=600, bbox_inches="tight", pad_inches=0.03,
                pil_kwargs={"compression": "tiff_lzw"})
    fig.savefig(FIG / f"{name}.pdf", bbox_inches="tight", pad_inches=0.03)
    fig.savefig(FIG / f"{name}.svg", bbox_inches="tight", pad_inches=0.03)


def read_csv(path):
    with path.open(encoding="utf-8-sig", newline="") as f:
        return list(csv.DictReader(f))


def holm_adjust(pvalues):
    pvalues = np.asarray(pvalues, float)
    order = np.argsort(pvalues)
    adjusted = np.empty_like(pvalues)
    running = 0.0
    m = len(pvalues)
    for rank, idx in enumerate(order):
        value = min(1.0, (m - rank) * pvalues[idx])
        running = max(running, value)
        adjusted[idx] = running
    return adjusted


def bootstrap_ci(values, n_boot=20000, seed=20261113):
    values = np.asarray(values, float)
    rng = np.random.default_rng(seed)
    means = np.empty(n_boot)
    for i in range(n_boot):
        means[i] = rng.choice(values, values.size, replace=True).mean()
    return np.percentile(means, [2.5, 97.5])


def rank_biserial(diff):
    diff = np.asarray(diff, float)
    diff = diff[np.abs(diff) > 1e-12]
    if diff.size == 0:
        return 0.0
    ranks = np.argsort(np.argsort(np.abs(diff))) + 1.0
    positive = ranks[diff > 0].sum()
    negative = ranks[diff < 0].sum()
    return float((positive - negative) / (positive + negative))


def load_combined():
    old = read_csv(OLD / "all_study_runs.csv")
    new = read_csv(RAW / "transfer_runs.csv")
    rows = []
    for r in old:
        variant = r["variant"]
        if r["study"] == "comparison" and variant == "CPO":
            label = "CPO + 2-opt"
        elif r["study"] == "comparison" and variant == "IPSO":
            label = "IPSO"
        elif r["study"] == "comparison" and variant == "CPO-DBO":
            label = "CPO-DBO + 2-opt"
        elif r["study"] == "ablation" and variant == "No 2-opt":
            label = "CPO-DBO"
        else:
            continue
        rows.append({
            "variant": label,
            "instance": int(r["instance"]),
            "repeat": int(r["repeat"]),
            "seed": int(r["seed"]),
            "objective": float(r["objective_au"]),
            "hash": r["initial_population_sha256"],
        })
    for r in new:
        label = "CPO" if r["variant"] == "CPO" else "IPSO + 2-opt"
        rows.append({
            "variant": label,
            "instance": int(r["instance"]),
            "repeat": int(r["repeat"]),
            "seed": int(r["seed"]),
            "objective": float(r["objective_au"]),
            "hash": r["initial_population_sha256"],
        })
    return rows


def main():
    apply_style()
    rows = load_combined()
    variants = ["CPO", "CPO + 2-opt", "IPSO", "IPSO + 2-opt", "CPO-DBO", "CPO-DBO + 2-opt"]
    values = {v: np.array([r["objective"] for r in rows if r["variant"] == v]) for v in variants}
    assert all(x.size == 200 for x in values.values())

    # Verify paired initialization hashes for the two newly added transfer controls.
    by_key = {(r["variant"], r["instance"], r["repeat"]): r for r in rows}
    for base, plus in (("CPO", "CPO + 2-opt"), ("IPSO", "IPSO + 2-opt")):
        for instance in range(1, 21):
            for repeat in range(1, 11):
                assert by_key[(base, instance, repeat)]["hash"] == by_key[(plus, instance, repeat)]["hash"]

    instance_means = {}
    for v in variants:
        instance_means[v] = np.array([
            np.mean([r["objective"] for r in rows if r["variant"] == v and r["instance"] == i])
            for i in range(1, 21)
        ])

    fried = friedmanchisquare(*(instance_means[v] for v in variants))

    pairs = [
        ("CPO", "CPO + 2-opt", "CPO"),
        ("IPSO", "IPSO + 2-opt", "IPSO"),
        ("CPO-DBO", "CPO-DBO + 2-opt", "CPO-DBO"),
    ]
    tests = []
    raw_p = []
    for no, yes, family in pairs:
        diff = instance_means[no] - instance_means[yes]
        w = wilcoxon(diff, alternative="two-sided", method="auto")
        ci = bootstrap_ci(diff)
        raw_p.append(float(w.pvalue))
        tests.append({
            "framework": family,
            "mean_without": float(values[no].mean()),
            "mean_with": float(values[yes].mean()),
            "mean_improvement": float(diff.mean()),
            "ci_low": float(ci[0]),
            "ci_high": float(ci[1]),
            "wilcoxon_statistic": float(w.statistic),
            "raw_p": float(w.pvalue),
            "rank_biserial": rank_biserial(diff),
            "wins": int(np.sum(diff > 0)),
            "ties": int(np.sum(np.abs(diff) <= 1e-12)),
            "losses": int(np.sum(diff < 0)),
        })
    adjusted = holm_adjust(raw_p)
    for row, p in zip(tests, adjusted):
        row["holm_p"] = float(p)

    # Difference-in-differences tests: is the 2-opt effect larger in CPO-DBO?
    effect_arrays = {
        "CPO": instance_means["CPO"] - instance_means["CPO + 2-opt"],
        "IPSO": instance_means["IPSO"] - instance_means["IPSO + 2-opt"],
        "CPO-DBO": instance_means["CPO-DBO"] - instance_means["CPO-DBO + 2-opt"],
    }
    interactions = []
    interaction_raw = []
    for comparator in ("CPO", "IPSO"):
        diff = effect_arrays["CPO-DBO"] - effect_arrays[comparator]
        w = wilcoxon(diff, alternative="two-sided", method="auto")
        ci = bootstrap_ci(diff, seed=30300 + len(interactions))
        interaction_raw.append(float(w.pvalue))
        interactions.append({
            "contrast": f"2-opt effect in CPO-DBO minus 2-opt effect in {comparator}",
            "mean_difference_in_effect": float(diff.mean()),
            "ci_low": float(ci[0]),
            "ci_high": float(ci[1]),
            "wilcoxon_p": float(w.pvalue),
            "rank_biserial": rank_biserial(diff),
        })
    for row, p in zip(interactions, holm_adjust(interaction_raw)):
        row["holm_p"] = float(p)

    # CPO-DBO reference tests across all five alternatives.
    reference = "CPO-DBO + 2-opt"
    refs = []
    ref_raw = []
    for comp in [v for v in variants if v != reference]:
        diff = instance_means[comp] - instance_means[reference]
        w = wilcoxon(diff, alternative="two-sided", method="auto")
        ci = bootstrap_ci(diff, seed=20261113 + len(refs))
        ref_raw.append(float(w.pvalue))
        refs.append({
            "reference": reference,
            "comparator": comp,
            "mean_difference_comparator_minus_reference": float(diff.mean()),
            "ci_low": float(ci[0]),
            "ci_high": float(ci[1]),
            "wilcoxon_p": float(w.pvalue),
            "rank_biserial": rank_biserial(diff),
        })
    ref_adj = holm_adjust(ref_raw)
    for row, p in zip(refs, ref_adj):
        row["holm_p"] = float(p)

    # Expanded nine-algorithm comparison requested by the reviewer.
    old_rows = read_csv(OLD / "all_study_runs.csv")
    expanded_names = [
        "CPO-DBO + 2-opt", "CPO", "CPO + 2-opt", "IPSO", "IPSO + 2-opt",
        "ACO-GA", "GA-PSO", "DBO", "GA",
    ]
    expanded_means = {name: instance_means[name] for name in
                      ("CPO-DBO + 2-opt", "CPO", "CPO + 2-opt", "IPSO", "IPSO + 2-opt")}
    expanded_run_values = {name: values[name] for name in expanded_means}
    for name in ("ACO-GA", "GA-PSO", "DBO", "GA"):
        selected = [r for r in old_rows if r["study"] == "comparison" and r["variant"] == name]
        expanded_run_values[name] = np.array([float(r["objective_au"]) for r in selected])
        expanded_means[name] = np.array([
            np.mean([float(r["objective_au"]) for r in selected if int(r["instance"]) == i])
            for i in range(1, 21)
        ])
    expanded_fried = friedmanchisquare(*(expanded_means[name] for name in expanded_names))
    rank_matrix = np.column_stack([expanded_means[name] for name in expanded_names])
    average_ranks = np.mean(np.apply_along_axis(rankdata, 1, rank_matrix), axis=0)
    expanded_tests = []
    expanded_raw = []
    for name in expanded_names[1:]:
        diff = expanded_means[name] - expanded_means[reference]
        w = wilcoxon(diff, alternative="two-sided", method="auto")
        ci = bootstrap_ci(diff, seed=41000 + len(expanded_tests))
        expanded_raw.append(float(w.pvalue))
        expanded_tests.append({
            "reference": reference,
            "comparator": name,
            "mean_difference_comparator_minus_reference": float(diff.mean()),
            "ci_low": float(ci[0]),
            "ci_high": float(ci[1]),
            "wilcoxon_p": float(w.pvalue),
            "rank_biserial": rank_biserial(diff),
        })
    for row, p in zip(expanded_tests, holm_adjust(expanded_raw)):
        row["holm_p"] = float(p)
    expanded_summary = []
    expanded_p = {r["comparator"]: r["holm_p"] for r in expanded_tests}
    for idx, name in enumerate(expanded_names):
        x = expanded_run_values[name]
        expanded_summary.append({
            "algorithm": name,
            "mean": float(x.mean()),
            "sd": float(x.std(ddof=1)),
            "median": float(np.median(x)),
            "instance_average_rank": float(average_ranks[idx]),
            "holm_p_vs_cpo_dbo": "" if name == reference else float(expanded_p[name]),
        })

    summary = [{
        "variant": v,
        "mean": float(values[v].mean()),
        "sd": float(values[v].std(ddof=1)),
        "median": float(np.median(values[v])),
        "instance_mean_ci_low": float(bootstrap_ci(instance_means[v], seed=20261113 + i)[0]),
        "instance_mean_ci_high": float(bootstrap_ci(instance_means[v], seed=20261113 + i)[1]),
    } for i, v in enumerate(variants)]

    def write_rows(name, data):
        with (ROOT / name).open("w", newline="", encoding="utf-8-sig") as f:
            w = csv.DictWriter(f, fieldnames=list(data[0]))
            w.writeheader(); w.writerows(data)

    write_rows("twoopt_transfer_summary.csv", summary)
    write_rows("twoopt_transfer_pairwise.csv", tests)
    write_rows("twoopt_transfer_interactions.csv", interactions)
    write_rows("expanded_reference_tests.csv", refs)
    write_rows("expanded_nine_algorithm_summary.csv", expanded_summary)
    write_rows("expanded_nine_algorithm_tests.csv", expanded_tests)
    (ROOT / "twoopt_transfer_statistics.json").write_text(json.dumps({
        "friedman_chi2": float(fried.statistic),
        "friedman_p": float(fried.pvalue),
        "summary": summary,
        "transfer_tests": tests,
        "interaction_tests": interactions,
        "reference_tests": refs,
        "expanded_nine_algorithm_friedman_chi2": float(expanded_fried.statistic),
        "expanded_nine_algorithm_friedman_p": float(expanded_fried.pvalue),
        "expanded_nine_algorithm_summary": expanded_summary,
        "expanded_nine_algorithm_tests": expanded_tests,
    }, indent=2), encoding="utf-8")

    # Publication figure: common legend above titles, no grid, four-sided frame.
    fig, axes = plt.subplots(1, 2, figsize=(7.2, 3.35))
    fig.subplots_adjust(left=0.09, right=0.985, bottom=0.19, top=0.69, wspace=0.34)
    ax = axes[0]
    frameworks = ["CPO", "IPSO", "CPO-DBO"]
    xpos = np.arange(3)
    means_no = [values[x].mean() for x in ("CPO", "IPSO", "CPO-DBO")]
    means_yes = [values[x].mean() for x in ("CPO + 2-opt", "IPSO + 2-opt", "CPO-DBO + 2-opt")]
    ci_no = [bootstrap_ci(instance_means[x], seed=101+i) for i, x in enumerate(("CPO", "IPSO", "CPO-DBO"))]
    ci_yes = [bootstrap_ci(instance_means[x], seed=201+i) for i, x in enumerate(("CPO + 2-opt", "IPSO + 2-opt", "CPO-DBO + 2-opt"))]
    for x, mn, my, cn, cy in zip(xpos, means_no, means_yes, ci_no, ci_yes):
        ax.plot([x-0.12, x+0.12], [mn, my], color="#A7A7A7", lw=1.0, zorder=1)
        ax.errorbar(x-0.12, mn, yerr=[[mn-cn[0]], [cn[1]-mn]], fmt="o", ms=5.5,
                    color=GRAY, mfc="white", mew=1.2, capsize=3, lw=1.1, zorder=3)
        ax.errorbar(x+0.12, my, yerr=[[my-cy[0]], [cy[1]-my]], fmt="s", ms=5.2,
                    color=ORANGE, mfc=ORANGE, mew=1.0, capsize=3, lw=1.1, zorder=3)
        label_y = max(cn[1], cy[1]) + 1.0
        ax.text(x, label_y, f"{mn-my:+.2f}", ha="center", va="bottom", fontsize=7.8, color="black")
    ax.set_xticks(xpos, frameworks)
    ax.set_ylabel("Relative exposure cost (a.u.)", labelpad=6)
    ax.set_title("Effect of transferring 2-opt to strong baselines", pad=8)
    ax.set_ylim(345, 380)
    style_axis(ax)

    ax = axes[1]
    effects = np.array([x["mean_improvement"] for x in tests])
    lows = np.array([x["ci_low"] for x in tests])
    highs = np.array([x["ci_high"] for x in tests])
    y = np.arange(3)
    ax.axvline(0, color="#777777", ls="--", lw=1.0, zorder=0)
    for i, (eff, lo, hi, row) in enumerate(zip(effects, lows, highs, tests)):
        color = GREEN if row["holm_p"] < 0.05 else BLUE
        ax.errorbar(eff, i, xerr=[[eff-lo], [hi-eff]], fmt="o", ms=5.8,
                    color=color, mfc=color, capsize=3, lw=1.3, zorder=3)
        ptxt = "<0.001" if row["holm_p"] < 0.001 else f"={row['holm_p']:.3f}"
        ax.text(hi + 0.5, i, f"$p_{{Holm}}${ptxt}", ha="left", va="center", fontsize=7.6)
    ax.set_yticks(y, frameworks)
    ax.invert_yaxis()
    ax.set_xlabel("Mean improvement from 2-opt (a.u.)", labelpad=6)
    ax.set_title("Paired instance-level effects (95% bootstrap CI)", pad=8)
    ax.set_xlim(min(-3.5, lows.min()-1), highs.max()+6.5)
    style_axis(ax)

    handles = [
        mpl.lines.Line2D([], [], color=GRAY, marker="o", mfc="white", lw=0, label="Without 2-opt"),
        mpl.lines.Line2D([], [], color=ORANGE, marker="s", mfc=ORANGE, lw=0, label="With 2-opt"),
    ]
    fig.legend(handles=handles, loc="upper center", bbox_to_anchor=(0.5, 0.985), ncol=2,
               frameon=False, handletextpad=0.6, columnspacing=2.2)
    fig.text(0.018, 0.705, "(a)", fontweight="bold", fontsize=10.5, ha="left", va="bottom")
    fig.text(0.515, 0.705, "(b)", fontweight="bold", fontsize=10.5, ha="left", va="bottom")
    save_all(fig, "Fig_transfer_2opt")
    plt.close(fig)

    # Corrected Fig. 21 based on the revised edge-level source of truth.
    src = Path(__file__).resolve().parents[1] / "data" / "obstacle_case3" / "path_edge_detail_mixed.csv"
    edge_rows = read_csv(src)
    totals = {"Dose": 0.0, "Length": 0.0}
    by_type = {"LOS": {"Dose": 0.0, "Length": 0.0}, "A*": {"Dose": 0.0, "Length": 0.0}}
    for r in edge_rows:
        typ = r["EdgeType"]
        for metric in ("Dose", "Length"):
            val = float(r[metric])
            totals[metric] += val
            by_type[typ][metric] += val
    los = np.array([100*by_type["LOS"][m]/totals[m] for m in ("Dose", "Length")])
    astar = 100 - los
    fig, ax = plt.subplots(figsize=(3.55, 3.05))
    fig.subplots_adjust(left=0.18, right=0.97, bottom=0.19, top=0.71)
    x = np.arange(2)
    ax.bar(x, los, width=0.52, color=BLUE, edgecolor="black", linewidth=0.6, label="LOS edges")
    ax.bar(x, astar, width=0.52, bottom=los, color=ORANGE, edgecolor="black", linewidth=0.6,
           label="A* detour edges")
    for i in range(2):
        ax.text(i, los[i]/2, f"{los[i]:.1f}%", ha="center", va="center", color="white", fontsize=8.5)
        ax.text(i, los[i]+astar[i]/2, f"{astar[i]:.1f}%", ha="center", va="center", color="black", fontsize=8.5)
    ax.set_xticks(x, ["Exposure cost", "Distance"])
    ax.set_ylabel("Share of tour total (%)", labelpad=6)
    ax.set_ylim(0, 105)
    ax.set_title("Edge composition of the obstacle-aware tour", pad=8)
    style_axis(ax)
    fig.legend(loc="upper center", bbox_to_anchor=(0.57, 0.985), ncol=2, frameon=False,
               handlelength=1.8, columnspacing=1.5)
    save_all(fig, "Fig21_edge_composition_corrected")
    plt.close(fig)

    # Corrected Fig. 22 assembled from the archived path panels. Repeated legends
    # and manually rounded titles are replaced by one shared legend and values
    # rounded from the common result source.
    source_image = Image.open(Path(__file__).resolve().parents[1] / "figures" / "legacy" / "image145.png").convert("RGB")
    crops = [
        source_image.crop((215, 44, 722, 526)),
        source_image.crop((915, 44, 1435, 526)),
        source_image.crop((215, 613, 722, 1095)),
        source_image.crop((915, 613, 1435, 1095)),
    ]
    titles = [
        "Cost = 216.78 a.u.; distance = 1070.1; A* edges = 2",
        "Cost = 215.57 a.u.; distance = 1060.7; A* edges = 3",
        "Cost = 216.78 a.u.; distance = 1064.7; A* edges = 3",
        "Cost = 215.41 a.u.; distance = 1057.8; A* edges = 2",
    ]
    fig, axes = plt.subplots(2, 2, figsize=(7.2, 6.3))
    fig.subplots_adjust(left=0.055, right=0.985, bottom=0.065, top=0.75,
                        wspace=0.07, hspace=0.30)
    for i, (ax, panel, title) in enumerate(zip(axes.ravel(), crops, titles)):
        ax.imshow(panel)
        ax.set_axis_off()
        ax.set_title(title, fontsize=8.5, pad=5, y=1.01)
        ax.text(-0.015, 1.13, f"({chr(97+i)})", transform=ax.transAxes,
                fontsize=10.5, fontweight="bold", ha="left", va="bottom")
    legend_handles = [
        mpl.lines.Line2D([], [], color="black", marker="s", mfc="#777777", lw=0,
                         ms=5, label="Obstacle"),
        mpl.lines.Line2D([], [], color="black", marker="o", mfc="black", lw=0,
                         ms=4, label="Inspection node"),
        mpl.lines.Line2D([], [], color=RED, marker="*", mfc=RED, lw=0,
                         ms=7, label="Radiation source"),
        mpl.lines.Line2D([], [], color="#D55E00", lw=1.4, label="LOS edge"),
        mpl.lines.Line2D([], [], color="#009E73", lw=1.4, ls="--", label="A* detour edge"),
    ]
    fig.legend(handles=legend_handles, loc="upper center", bbox_to_anchor=(0.5, 0.985),
               ncol=3, frameon=False, handlelength=2.0, columnspacing=1.6)
    fig.text(0.52, 0.012, "X coordinate (m)", ha="center", va="bottom", fontsize=9.5)
    fig.text(0.012, 0.42, "Y coordinate (m)", ha="left", va="center",
             rotation=90, fontsize=9.5)
    save_all(fig, "Fig22_obstacle_combinations_corrected")
    plt.close(fig)

    print(json.dumps({"friedman": [fried.statistic, fried.pvalue], "tests": tests}, indent=2))


if __name__ == "__main__":
    main()
