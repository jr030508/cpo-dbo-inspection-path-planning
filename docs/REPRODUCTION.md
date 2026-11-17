# Reproduction guide

All commands are executed from the `python/` directory.
Every experiment writes to `output/<script_name>/`; the committed copies in `results/`
are the authoritative outputs referenced by the manuscript.

## 0. Environment

```bash
pip install -r requirements.txt   # numpy==2.5.4, numba==0.68.1, scipy, pandas, matplotlib, Pillow
```

Reference machine: Intel Core i7-14650HX, 32 GB RAM, Windows 11 64-bit,
Python 3.12.14. Experiments run single-process, single-thread
(`OMP_NUM_THREADS` etc. are forced to 1 inside the scripts).

## 1. Instance library

The 20 perturbed 51-node instances are committed under `data/instance_lib20/`
(one `doseMatrix_instXX.csv` per instance). They were generated deterministically
by MATLAB:

```matlab
cd data/instance_generator
build_instance_library20_case3   % rng(6,'twister'), reproduces identical matrices
```

Perturbation rules (identical for every instance): source strength ±10%
(renormalized), source positions ±2 m, 1–3% symmetric noise on the dose matrix,
top-5% dose edges amplified ×3.

## 2. Validation

```bash
python validate_suite.py
```

Loads `data/instance_lib20/doseMatrix_inst01.csv` and sanity-checks the solver
suite (budget accounting, permutation feasibility).

## 3. Nine-configuration comparison, ablation, sensitivity

```bash
python run_studies.py        # heavy: 9 configurations × 20 instances × 10 seeds
python analyze_studies.py    # recomputes objectives from raw runs, all statistical tests
```

Protocol: 1,960,400 complete route evaluations per run; paired seeds
(seed = 20000 + 100*instance + repeat); identical initial populations within pairs.
`analyze_studies.py` asserts 4,000 runs and recomputation error < 1e-9, then emits
Friedman tests, Holm-corrected Wilcoxon tests, instance-level bootstrap CIs,
rank-biserial effect sizes.

Outputs (also committed under `results/supplementary/` and `results/raw_runs/`):
`study_summary.csv`, `all_study_runs.csv`, `statistics.json`,
`algorithm_pairwise_tests.csv`, `ablation_pairwise_tests.csv`,
`sensitivity_*.csv`, `instance_level_means.csv`.

Key manuscript numbers reproduced here: CPO-DBO mean 359.86 a.u., mean instance
rank 3.00, χ²(8) = 64.147 (p < 0.001); ablation deltas 6.44 / 6.15 a.u.
(Holm p = 0.027 / 0.040).

## 4. 2-opt transfer experiment

```bash
python run_transfer_experiment.py     # CPO, IPSO, CPO-DBO(−2-opt) with sampled 2-opt
python analyze_and_plot_transfer.py   # pairwise CIs, effect-difference tests, Fig. 18
```

The third row of the committed `results/transfer/twoopt_transfer_pairwise.csv`
(6.44 a.u., CI [2.73, 10.21], Holm p = 0.017) is the transfer-experiment bootstrap
of the same paired data as the ablation test (CI [2.70, 10.17] in the ablation
table comes from the ablation-family bootstrap run).

## 5. Optimizer-only timing

```bash
python run_single_process_timing.py
```

Measures the optimizer core only: one explicit JIT warm-up per solver, then 200
timed optimizer calls in randomized order, single process and thread. Excludes
dose-matrix construction, CSV I/O, serialization, and plotting. Median results:
0.779 s (2-opt kept) vs. 0.694 s (2-opt removed).

## 6. Legacy MATLAB figures

```bash
python redraw_legacy_figures.py   # reads data/benchmarks/*/results_*.mat
```

Regenerates the descriptive MATLAB-era figures (box plots, ECDF, convergence) in
publication style from the committed benchmark outputs.

## 6b. Authoritative paper-figure sources

The final manuscript figures (labels: *relative exposure cost (a.u.*) are in
`figures/paper/`, and their editable MATLAB `.fig` sources in
`figures/fig_sources/` — see `figures/README.md` for the figure-number map.

Three residual label issues are tracked in `matlab/figure_label_fixes/`:
run `fix_labels` there to (1) rename the 图6 legend "Optimal Dose" →
"Optimal exposure cost", (2) fix the 图4 Y-axis typo "Start Npde" →
"Start Node", and (3) print the 图21(a) title so its value (212.364 a.u.)
can be reconciled with Table 15 (212.35 a.u.) by hand.

## 7. Standardized obstacle re-run (mechanism illustration only)

```bash
python run_obstacle_standardized.py
```

Single fixed seed (6) — the nominal case used for the Section 3.7 obstacle
mechanism analysis. This is a mechanism demonstration, not a statistical result;
see the manuscript for the scope statement.
