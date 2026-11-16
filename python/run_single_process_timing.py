import csv
import ctypes
import json
import os
import platform
import random
import sys
import time
from pathlib import Path

for key in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS", "NUMBA_NUM_THREADS"):
    os.environ[key] = "1"

import numba
import numpy as np

BASE = Path(__file__).resolve().parent
sys.path.insert(0, str(BASE))
from solver_suite import BUDGET, run_cpo, run_cpo_dbo, run_ipso  # noqa: E402

PYDIR = Path(__file__).resolve().parent
OUT = Path(__file__).resolve().parents[1] / "output" / "run_single_process_timing"
OUT.mkdir(parents=True, exist_ok=True)
sys.path.insert(0, str(PYDIR))
from run_transfer_experiment import run_cpo_no2opt, run_ipso_plus2opt  # noqa: E402

DATA = Path(__file__).resolve().parents[1] / "data" / "instance_lib20"


def memory_gb():
    class MEMORYSTATUSEX(ctypes.Structure):
        _fields_ = [
            ("dwLength", ctypes.c_ulong),
            ("dwMemoryLoad", ctypes.c_ulong),
            ("ullTotalPhys", ctypes.c_ulonglong),
            ("ullAvailPhys", ctypes.c_ulonglong),
            ("ullTotalPageFile", ctypes.c_ulonglong),
            ("ullAvailPageFile", ctypes.c_ulonglong),
            ("ullTotalVirtual", ctypes.c_ulonglong),
            ("ullAvailVirtual", ctypes.c_ulonglong),
            ("ullAvailExtendedVirtual", ctypes.c_ulonglong),
        ]
    status = MEMORYSTATUSEX()
    status.dwLength = ctypes.sizeof(status)
    ctypes.windll.kernel32.GlobalMemoryStatusEx(ctypes.byref(status))
    return status.ullTotalPhys / (1024**3)


def cpu_name():
    try:
        import winreg
        key = winreg.OpenKey(winreg.HKEY_LOCAL_MACHINE,
                             r"HARDWARE\DESCRIPTION\System\CentralProcessor\0")
        return winreg.QueryValueEx(key, "ProcessorNameString")[0].strip()
    except Exception:
        return platform.processor()


def execute(label, dose, seed):
    if label == "CPO":
        return run_cpo_no2opt(dose, seed)
    if label == "CPO + 2-opt":
        return run_cpo(dose, seed)
    if label == "IPSO":
        return run_ipso(dose, seed)
    if label == "IPSO + 2-opt":
        return run_ipso_plus2opt(dose, seed)
    if label == "CPO-DBO":
        return run_cpo_dbo(dose, seed, twoopt_attempts=0)
    if label == "CPO-DBO + 2-opt":
        return run_cpo_dbo(dose, seed, twoopt_attempts=1)
    raise ValueError(label)


def main():
    labels = ["CPO", "CPO + 2-opt", "IPSO", "IPSO + 2-opt", "CPO-DBO", "CPO-DBO + 2-opt"]
    load_times = {}
    matrices = {}
    for instance in range(1, 21):
        start = time.perf_counter()
        matrices[instance] = np.loadtxt(DATA / f"doseMatrix_inst{instance:02d}.csv", delimiter=",")
        load_times[instance] = time.perf_counter() - start

    # Explicit warm-up. These calls are not included in the reported timings.
    warm = matrices[1]
    for label in labels:
        execute(label, warm, 999001)
        print(f"warmed {label}", flush=True)

    schedule = [(instance, repeat, label)
                for instance in range(1, 21)
                for repeat in range(1, 11)
                for label in labels]
    random.Random(20261113).shuffle(schedule)
    rows = []
    started = time.perf_counter()
    for index, (instance, repeat, label) in enumerate(schedule, 1):
        seed = 20000 + 100 * instance + repeat
        start = time.perf_counter()
        result = execute(label, matrices[instance], seed)
        elapsed = time.perf_counter() - start
        rows.append({
            "variant": label,
            "instance": instance,
            "repeat": repeat,
            "seed": seed,
            "optimizer_only_seconds": elapsed,
            "objective_au": float(result[0]),
            "matrix_load_seconds": load_times[instance],
            "full_tour_evaluations": int(result[4]),
        })
        if index % 50 == 0:
            print(f"{index}/{len(schedule)} timed runs; elapsed={time.perf_counter()-started:.1f}s", flush=True)

    rows.sort(key=lambda x: (x["variant"], x["instance"], x["repeat"]))
    with (OUT / "single_process_timing_runs.csv").open("w", newline="", encoding="utf-8-sig") as f:
        writer = csv.DictWriter(f, fieldnames=list(rows[0]))
        writer.writeheader(); writer.writerows(rows)

    summary = []
    for label in labels:
        x = np.array([r["optimizer_only_seconds"] for r in rows if r["variant"] == label])
        summary.append({
            "variant": label,
            "n": int(x.size),
            "mean_seconds": float(x.mean()),
            "sd_seconds": float(x.std(ddof=1)),
            "median_seconds": float(np.median(x)),
            "q1_seconds": float(np.percentile(x, 25)),
            "q3_seconds": float(np.percentile(x, 75)),
            "min_seconds": float(x.min()),
            "max_seconds": float(x.max()),
        })
    with (OUT / "single_process_timing_summary.csv").open("w", newline="", encoding="utf-8-sig") as f:
        writer = csv.DictWriter(f, fieldnames=list(summary[0]))
        writer.writeheader(); writer.writerows(summary)

    manifest = {
        "scope": "optimizer-only time after explicit JIT warm-up",
        "excluded": ["dose-matrix construction", "CSV loading", "result serialization", "plotting", "JIT warm-up"],
        "processes": 1,
        "threads": 1,
        "run_order": "deterministically randomized across algorithms",
        "instances": 20,
        "seeds_per_instance": 10,
        "runs_per_algorithm": 200,
        "budget": BUDGET,
        "cpu": cpu_name(),
        "ram_gb": memory_gb(),
        "os": platform.platform(),
        "python": sys.version,
        "numpy": np.__version__,
        "numba": numba.__version__,
    }
    (OUT / "single_process_timing_manifest.json").write_text(
        json.dumps(manifest, indent=2), encoding="utf-8"
    )
    print(json.dumps({"summary": summary, "manifest": manifest}, indent=2), flush=True)


if __name__ == "__main__":
    main()
