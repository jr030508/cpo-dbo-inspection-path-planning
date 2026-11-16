import csv
import hashlib
import json
import os
import platform
import sys
import time
from concurrent.futures import ProcessPoolExecutor, as_completed
from pathlib import Path

for key in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS"):
    os.environ[key] = "1"

import numpy as np

from solver_suite import (
    BUDGET, fitness, run_aco_ga, run_cpo, run_cpo_dbo, run_dbo,
    run_ga, run_ga_pso, run_ipso,
)


OUT = Path(__file__).resolve().parents[1] / "output" / Path(__file__).stem
OUT.mkdir(parents=True, exist_ok=True)
ROOT = Path(__file__).resolve().parents[1] / "data"
RUNS = OUT / "runs"
RUNS.mkdir(exist_ok=True)

CONFIGS = [
    ("CPO-DBO", "comparison", "cpo_dbo", {}),
    ("CPO", "comparison", "cpo", {}),
    ("DBO", "comparison", "dbo", {}),
    ("GA", "comparison", "ga", {}),
    ("IPSO", "comparison", "ipso", {}),
    ("GA-PSO", "comparison", "ga_pso", {}),
    ("ACO-GA", "comparison", "aco_ga", {}),
    ("No 2-opt", "ablation", "cpo_dbo", {"twoopt_attempts": 0}),
    ("Fixed probabilities", "ablation", "cpo_dbo", {"dynamic_probs": 0}),
    ("No group bias", "ablation", "cpo_dbo", {"use_group_bias": 0}),
    ("No elitism", "ablation", "cpo_dbo", {"use_elite": 0}),
    ("Swap-only mutation", "ablation", "cpo_dbo", {"mutation_mode": 1}),
    ("No personal-best guidance", "ablation", "cpo_dbo", {"use_pbest": 0}),
    ("M=200", "sensitivity", "cpo_dbo", {"m": 200}),
    ("M=600", "sensitivity", "cpo_dbo", {"m": 600}),
    ("Tscale=5000", "sensitivity", "cpo_dbo", {"tscale": 5000.0}),
    ("Tscale=15000", "sensitivity", "cpo_dbo", {"tscale": 15000.0}),
    ("Elite=1%", "sensitivity", "cpo_dbo", {"elite_ratio": 0.01}),
    ("Elite=5%", "sensitivity", "cpo_dbo", {"elite_ratio": 0.05}),
    ("2-opt=3", "sensitivity", "cpo_dbo", {"twoopt_attempts": 3}),
]


def safe_label(label):
    return (label.replace(" ", "_").replace("=", "-").replace("%", "pct")
            .replace("/", "_").replace("–", "-").replace("-", "_"))


def execute(kind, dose, seed, kwargs):
    if kind == "cpo_dbo":
        return run_cpo_dbo(dose, seed, budget=BUDGET, **kwargs)
    if kind == "cpo":
        return run_cpo(dose, seed, budget=BUDGET)
    if kind == "dbo":
        return run_dbo(dose, seed, budget=BUDGET)
    if kind == "ga":
        return run_ga(dose, seed, budget=BUDGET)
    if kind == "ipso":
        return run_ipso(dose, seed, budget=BUDGET)
    if kind == "ga_pso":
        return run_ga_pso(dose, seed, budget=BUDGET)
    if kind == "aco_ga":
        return run_aco_ga(dose, seed, budget=BUDGET)
    raise ValueError(kind)


def run_pair(instance, repeat):
    dose = np.loadtxt(ROOT / "instance_lib20" / f"doseMatrix_inst{instance:02d}.csv", delimiter=",")
    seed = 20000 + 100 * instance + repeat
    records = []
    for label, study, kind, kwargs in CONFIGS:
        folder = RUNS / safe_label(label)
        folder.mkdir(exist_ok=True)
        destination = folder / f"inst{instance:02d}_seed{seed}.npz"
        if destination.exists():
            with np.load(destination) as saved:
                records.append(json.loads(str(saved["record"])))
            continue
        start = time.perf_counter()
        result = execute(kind, dose, seed, kwargs)
        elapsed = time.perf_counter() - start
        value, path, curve, first_hit, evaluations = result[:5]
        initial = result[-1]
        full = np.r_[0, path, 0]
        recomputed = float(sum(dose[a, b] for a, b in zip(full[:-1], full[1:])))
        assert np.array_equal(np.sort(path), np.arange(1, 51))
        assert abs(recomputed - value) < 1e-9
        assert np.all(np.diff(curve) <= 1e-10)
        record = {
            "study": study,
            "variant": label,
            "instance": instance,
            "repeat": repeat,
            "seed": seed,
            "objective_au": float(value),
            "seconds": elapsed,
            "first_hit_iteration": int(first_hit),
            "full_tour_evaluations": int(evaluations),
            "initial_population_sha256": hashlib.sha256(initial.tobytes()).hexdigest(),
            "recomputed_objective_au": recomputed,
        }
        if kind == "cpo_dbo":
            record["two_opt_attempts"] = int(result[5])
            record["two_opt_acceptances"] = int(result[6])
        else:
            record["two_opt_attempts"] = 0
            record["two_opt_acceptances"] = 0
        np.savez_compressed(
            destination,
            path=np.r_[1, path + 1, 1],
            convergence_progress=curve,
            record=json.dumps(record),
        )
        records.append(record)
    return records


def main():
    manifest = {
        "engine": "Python equivalent permutation implementation; not bitwise MATLAB reproduction",
        "objective_unit": "relative exposure cost (a.u.)",
        "budget": BUDGET,
        "instances": 20,
        "paired_seeds_per_instance": 10,
        "python": sys.version,
        "numpy": np.__version__,
        "platform": platform.platform(),
        "workers": 4,
        "seed_formula": "20000 + 100*instance + repeat",
        "configs": [{"label": a, "study": b, "kind": c, "kwargs": d} for a, b, c, d in CONFIGS],
    }
    (OUT / "experiment_manifest.json").write_text(
        json.dumps(manifest, ensure_ascii=False, indent=2), encoding="utf-8"
    )
    rows = []
    started = time.perf_counter()
    with ProcessPoolExecutor(max_workers=4) as pool:
        jobs = [pool.submit(run_pair, instance, repeat)
                for instance in range(1, 21) for repeat in range(1, 11)]
        for completed, job in enumerate(as_completed(jobs), 1):
            rows.extend(job.result())
            (OUT / "progress.json").write_text(json.dumps({
                "completed_pairs": completed,
                "total_pairs": 200,
                "elapsed_seconds": time.perf_counter() - started,
            }), encoding="utf-8")
            if completed % 10 == 0:
                print(completed, "/ 200 pairs", round(time.perf_counter() - started, 1), "s", flush=True)
    rows.sort(key=lambda x: (x["study"], x["variant"], x["instance"], x["repeat"]))
    with (OUT / "all_study_runs.csv").open("w", newline="", encoding="utf-8-sig") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)
    print("COMPLETED", len(rows), "runs in", round(time.perf_counter() - started, 1), "s")


if __name__ == "__main__":
    main()
