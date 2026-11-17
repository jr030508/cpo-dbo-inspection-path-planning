# CPO-DBO Hybrid Metaheuristic for Discrete-Node Inspection Path Planning in Radioactive Scenarios

Code, data, and run-level results for the manuscript *"Radioactive-scenario discrete-node
inspection path planning based on a CPO-DBO hybrid metaheuristic"*.

The repository contains everything needed to reproduce the statistical experiments in the
paper: the dose-field / instance generation, the discrete solver suite (nine algorithm
configurations), the paired multi-instance experiment protocols, the 2-opt transfer
(gauge) experiment, single-process timing, and all statistical analyses
(Friedman + Holm-corrected Wilcoxon, exact sign-flip tests, bootstrap CIs,
rank-biserial effect sizes).

## Repository layout

```
├── data/
│   ├── instance_lib20/          20 perturbed 51-node dose matrices (rng(6), deterministic)
│   ├── instance_generator/      MATLAB generator for the instance library
│   ├── benchmarks/              MATLAB benchmark outputs (per-algorithm, 20 instances)
│   └── obstacle_case3/          LOS / mixed cost-dose-edgeType matrices, obstacle geometry
├── matlab/                      Legacy MATLAB implementation (descriptive experiments)
│   ├── single/                  IPSO / GA / DBO / CPO
│   ├── hybrid/                  CPO-DBO (30- and 51-node variants)
│   ├── benchmark/               20-instance benchmark scripts
│   └── obstacle/                LOS vs. mixed model analysis (Section 3.7)
├── python/                      Python/Numba implementation (all statistical experiments)
├── results/                     Committed outputs of every experiment
│   ├── supplementary/           Summary + test CSVs/JSON (seven-algorithm study, ablation, sensitivity)
│   ├── transfer/                2-opt transfer experiment (nine configurations)
│   ├── timing/                  Single-process timing
│   ├── obstacle/                Standardized obstacle re-run (mechanism illustration)
│   └── raw_runs/                Run-level records (4,000 runs + transfer + timing)
├── figures/                     Paper figures
│   ├── paper/                   Final exported figures (exposure-cost labels)
│   ├── fig_sources/             MATLAB .fig sources for every redrawn panel
│   ├── revision|supplementary/  Python-pipeline analysis figures
│   └── legacy/                  Earlier drafts (provenance only)
└── docs/REPRODUCTION.md         Step-by-step reproduction guide
```

## Environment

All Python experiments were executed on a single machine:

| Item | Value |
|---|---|
| CPU | Intel Core i7-14650HX |
| RAM | 32 GB |
| OS | Windows 11 64-bit |
| Python | 3.12.14 |
| NumPy | 2.5.4 |
| Numba | 0.68.1 |
| Execution mode | 1 process, 1 thread, randomized run order |

```bash
pip install -r requirements.txt
```

Legacy MATLAB experiments were run in MATLAB (see the manuscript for the descriptive
figures they produced; those results are descriptive and are not compared quantitatively
with the Python timings).

## Unified evaluation budget

Every algorithm configuration runs with the same budget of **1,960,400 complete
route evaluations per run** (population M = 400, ~4,900 generations), on the same
20-instance library, with 10 paired seeds per instance (200 runs per configuration).
Seed formula: `20000 + 100*instance + repeat`.

## Quick start

```bash
cd python

# 1) Validate the solver suite against the committed instance library
python validate_suite.py

# 2) Nine-configuration comparison + ablation + sensitivity (heavy: ~full protocol)
python run_studies.py
python analyze_studies.py

# 3) 2-opt transfer experiment (CPO / IPSO / CPO-DBO with and without 2-opt)
python run_transfer_experiment.py
python analyze_and_plot_transfer.py

# 4) Optimizer-only timing (JIT warm-up, single process/thread)
python run_single_process_timing.py

# 5) Redraw the legacy MATLAB descriptive figures from committed .mat outputs
python redraw_legacy_figures.py
```

Outputs are written to `output/` (git-ignored). The committed files under `results/`
are exactly the ones reported in the manuscript.

## Mapping to the manuscript

| Manuscript item | Source |
|---|---|
| Table 8 (nine-configuration statistics, χ²(8)=64.147) | `run_studies.py` + `expanded_nine_algorithm_*` / `analyze_studies.py` |
| Table 9–12 (ablation + sensitivity) | `run_studies.py` (variants) + `results/supplementary` |
| Table 13 / Fig. 18 (2-opt transfer) | `run_transfer_experiment.py` + `analyze_and_plot_transfer.py` |
| Timing rows (0.779 / 0.694 s medians) | `run_single_process_timing.py` |
| Instance library (Table 2 perturbation rules) | `data/instance_generator/build_instance_library20_case3.m` |
| Obstacle mechanism (Section 3.7) | `matlab/obstacle/` + `data/obstacle_case3/` |

## Notes on scope

* The standardized obstacle re-run (`run_obstacle_standardized.py`) is a single-seed
  mechanism illustration only and is **not** used for statistical claims in the manuscript.
* Objective values are relative transit-exposure costs in arbitrary units (a.u.),
  not calibrated physical doses.

## License

MIT — see [LICENSE](LICENSE).

## Citation

If you use this code, please cite the accompanying manuscript
(R. Jiang, H. Guan, Q. You, B. Cao — citation details to be updated upon acceptance).
