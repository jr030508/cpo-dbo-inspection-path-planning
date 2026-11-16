import csv
import hashlib
import json
import os
import sys
import time
from concurrent.futures import ProcessPoolExecutor, as_completed
from pathlib import Path

for key in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS"):
    os.environ[key] = "1"

import numpy as np
from numba import njit

BASE = Path(__file__).resolve().parent
sys.path.insert(0, str(BASE))
from solver_suite import (  # noqa: E402
    BUDGET,
    curve_checkpoints,
    fitness,
    initialize_population,
    evaluate_population,
    inject_segment,
    insert,
    invert,
    mutate,
    opt2,
    ox,
    scramble,
    swap,
)

ROOT = Path(__file__).resolve().parents[1] / "data"
OUT = Path(__file__).resolve().parents[1] / "output" / Path(__file__).stem
OUT.mkdir(parents=True, exist_ok=True)
RUNS = OUT / "runs_transfer"
RUNS.mkdir(parents=True, exist_ok=True)


@njit(cache=True)
def run_cpo_no2opt(dose, seed, m=400, budget=BUDGET):
    np.random.seed(seed)
    n = dose.shape[0] - 1
    iterations = int((budget - m) // m)
    pop, fits = initialize_population(m, n)
    evaluate_population(pop, dose, fits)
    initial = pop.copy()
    pbest = pop.copy()
    pfit = fits.copy()
    idx = np.argmin(fits)
    best = fits[idx]
    best_path = pop[idx].copy()
    curve = np.empty(iterations)
    for it in range(iterations):
        t = (it + 1) / iterations
        pd = 0.70 * (1.0 - t)
        pa = 0.30 * (1.0 - t)
        pe = 0.20 + 0.30 * t
        total = pd + pa + pe
        pd /= total
        pa /= total
        for i in range(m):
            r = np.random.random()
            if r < pd:
                child = ox(pop[i], best_path) if np.random.random() < 0.70 else pop[i].copy()
                if np.random.random() < 0.50:
                    swap(child)
                mutate(child, 0.60, 0.15, 0.15, 0.10)
            elif r < pd + pa:
                child = ox(pop[i], pbest[np.random.randint(m)])
                mutate(child, 0.35, 0.20, 0.25, 0.10)
            else:
                child = pop[i].copy()
                if np.random.random() < 0.85:
                    scramble(child)
                if np.random.random() < 0.65:
                    invert(child)
                mutate(child, 0.25, 0.25, 0.35, 0.25)
            value = fitness(child, dose)
            if value < fits[i]:
                pop[i] = child
                fits[i] = value
                if value < pfit[i]:
                    pbest[i] = child
                    pfit[i] = value
        idx = np.argmin(pfit)
        if pfit[idx] < best:
            best = pfit[idx]
            best_path = pbest[idx].copy()
        curve[it] = best
    first_hit = int(np.where(curve <= best + 1e-12)[0][0] + 1)
    return best, best_path, curve_checkpoints(curve), first_hit, m + iterations * m, initial


@njit(cache=True)
def run_ipso_plus2opt(dose, seed, m=400, budget=BUDGET):
    np.random.seed(seed)
    n = dose.shape[0] - 1
    iterations = int((budget - m) // m)
    chaos_iters = min(200, iterations)
    chaos = np.random.random((m, n))
    pop = np.empty((m, n), np.int64)
    fits = np.empty(m)
    attempts = 0
    accepted = 0
    for i in range(m):
        pop[i] = np.argsort(chaos[i]) + 1
        fits[i] = fitness(pop[i], dose)
    initial = pop.copy()
    pbest = pop.copy()
    pfit = fits.copy()
    idx = np.argmin(fits)
    best = fits[idx]
    best_path = pop[idx].copy()
    curve = np.empty(iterations)
    for it in range(chaos_iters):
        chaos = 4.0 * chaos * (1.0 - chaos)
        for i in range(m):
            child = np.argsort(chaos[i]) + 1
            attempts += 1
            accepted += opt2(child, dose, 1)
            value = fitness(child, dose)
            if value < pfit[i]:
                pbest[i] = child
                pfit[i] = value
        idx = np.argmin(pfit)
        if pfit[idx] < best:
            best = pfit[idx]
            best_path = pbest[idx].copy()
        curve[it] = best
    for it in range(chaos_iters, iterations):
        for i in range(m):
            child = ox(pop[i], pbest[i])
            child = ox(child, best_path)
            if np.random.random() < 0.10:
                child = inject_segment(child, best_path)
            mutate(child, 0.60, 0.15, 0.15, 0.10)
            attempts += 1
            accepted += opt2(child, dose, 1)
            value = fitness(child, dose)
            pop[i] = child
            fits[i] = value
            if value < pfit[i]:
                pbest[i] = child
                pfit[i] = value
        idx = np.argmin(pfit)
        if pfit[idx] < best:
            best = pfit[idx]
            best_path = pbest[idx].copy()
        curve[it] = best
    first_hit = int(np.where(curve <= best + 1e-12)[0][0] + 1)
    return (
        best,
        best_path,
        curve_checkpoints(curve),
        first_hit,
        m + iterations * m,
        attempts,
        accepted,
        initial,
    )


CONFIGS = [
    ("CPO", run_cpo_no2opt),
    ("IPSO+2-opt", run_ipso_plus2opt),
]


def safe_label(label):
    return label.replace("+", "_plus_").replace("-", "_").replace(" ", "_")


def run_pair(instance, repeat):
    dose = np.loadtxt(ROOT / "instance_lib20" / f"doseMatrix_inst{instance:02d}.csv", delimiter=",")
    seed = 20000 + 100 * instance + repeat
    records = []
    for label, solver in CONFIGS:
        folder = RUNS / safe_label(label)
        folder.mkdir(exist_ok=True)
        destination = folder / f"inst{instance:02d}_seed{seed}.npz"
        if destination.exists():
            with np.load(destination) as saved:
                records.append(json.loads(str(saved["record"])))
            continue
        start = time.perf_counter()
        result = solver(dose, seed)
        elapsed = time.perf_counter() - start
        value, path, curve, first_hit, evaluations = result[:5]
        initial = result[-1]
        full = np.r_[0, path, 0]
        recomputed = float(sum(dose[a, b] for a, b in zip(full[:-1], full[1:])))
        assert np.array_equal(np.sort(path), np.arange(1, 51))
        assert abs(recomputed - value) < 1e-9
        assert np.all(np.diff(curve) <= 1e-10)
        record = {
            "variant": label,
            "instance": instance,
            "repeat": repeat,
            "seed": seed,
            "objective_au": float(value),
            "seconds_parallel_run": elapsed,
            "first_hit_iteration": int(first_hit),
            "full_tour_evaluations": int(evaluations),
            "initial_population_sha256": hashlib.sha256(initial.tobytes()).hexdigest(),
            "recomputed_objective_au": recomputed,
            "two_opt_attempts": int(result[5]) if label == "IPSO+2-opt" else 0,
            "two_opt_acceptances": int(result[6]) if label == "IPSO+2-opt" else 0,
        }
        np.savez_compressed(
            destination,
            path=np.r_[1, path + 1, 1],
            convergence_progress=curve,
            record=json.dumps(record),
        )
        records.append(record)
    return records


def main():
    rows = []
    started = time.perf_counter()
    with ProcessPoolExecutor(max_workers=4) as pool:
        jobs = [pool.submit(run_pair, instance, repeat)
                for instance in range(1, 21) for repeat in range(1, 11)]
        for completed, job in enumerate(as_completed(jobs), 1):
            rows.extend(job.result())
            if completed % 10 == 0:
                print(f"{completed}/200 pairs; elapsed={time.perf_counter()-started:.1f}s", flush=True)
    rows.sort(key=lambda x: (x["variant"], x["instance"], x["repeat"]))
    with (OUT / "transfer_runs.csv").open("w", newline="", encoding="utf-8-sig") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)
    manifest = {
        "engine": "Python/Numba discrete permutation implementation",
        "new_variants": [x[0] for x in CONFIGS],
        "budget": BUDGET,
        "instances": 20,
        "paired_seeds_per_instance": 10,
        "workers_for_objective_experiment": 4,
        "seed_formula": "20000 + 100*instance + repeat",
        "python": sys.version,
        "numpy": np.__version__,
    }
    (OUT / "transfer_manifest.json").write_text(
        json.dumps(manifest, ensure_ascii=False, indent=2), encoding="utf-8"
    )
    print(f"COMPLETED {len(rows)} runs in {time.perf_counter()-started:.1f}s", flush=True)


if __name__ == "__main__":
    main()
