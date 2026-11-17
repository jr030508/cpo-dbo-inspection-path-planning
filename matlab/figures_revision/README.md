# MATLAB Figure Revision Scripts (relative exposure cost, a.u.)

The 7 figures in the original submission labelled optimal quantities with
"μSv" units. The model objective is a **relative transit exposure cost**
(integral of a dimensionless radiation field along the tour), so the labels
have been corrected to "Optimal exposure cost (a.u.)" with Times New Roman
typography. **No experiment was re-run**: all figure data come from the
archived result files in `data/` (identical to the ones used for the
original submission; values verified against Tables 7/15).

| Script | Output figure | Data source |
|---|---|---|
| `fig4_heatmap_cost.m`  | Fig. 4  30-node pairwise cost heatmap | original `.fig` / `fig4_dose30.mat` |
| `fig6_convergence_cost.m`  | Fig. 6  convergence curves (DBO/IPSO/CPO-DBO/CPO) | pixel-extracted curves, plateau anchored at 51.4888 |
| `fig8_boxplots_cost.m`  | Fig. 8  20-instance boxplots (a-d) | `data/bench_results/results_*_20instances.mat` |
| `fig10_ecdf_cost.m`  | Fig. 10 ECDF | same bench results |
| `fig11_risk_bars_cost.m`  | Fig. 11 P90/P95/Max | same bench results |
| `fig14_seed_curves_cost.m`  | Fig. 14 five seed curves + enlargement | pixel-extracted, plateau anchored at 327.719 |
| `fig21_free_vs_obstacle_cost.m`  | Fig. 21 free vs. obstacle best tours | `fig21_LOS.mat` / `fig21_MIXED.mat` (pathSum reproduces Table 15: 212.35/1046.5 and 216.78/1070.1) |

Run `RUN_ALL.m` to regenerate all 7 PNGs (300 dpi) into `output/`.

Figs 5, 9, 12, 13, 15-20, 22 and 23 already used correct a.u. labels and
are unaffected.
