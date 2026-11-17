# Paper figures

Two kinds of assets live here:

- `paper/` — final exported figures used in the manuscript (labels use
  *relative exposure cost (a.u.)*, Times New Roman).
- `fig_sources/` — the MATLAB `.fig` sources behind every redrawn panel, so
  any label/annotation can be adjusted and re-exported (`openfig` in MATLAB).

Legacy figures (earlier drafts, including dose-labelled versions that were
renamed to relative exposure cost) are kept in `figures/legacy/` for
provenance only — **do not use them in the manuscript**.

## Figure map (paper figure number -> file)

| Paper | File(s) in `paper/` | Source `.fig` |
|---|---|---|
| 图3  30-node base scene | (legacy `figures/legacy/`, unchanged) | — |
| 图4  Exposure-cost matrix heatmap (30 nodes) | `fig04_exposure_cost_heatmap_30nodes.jpg` | `fig04_exposure_cost_heatmap.fig` |
| 图6  Convergence curves (a)DBO (b)IPSO (c)CPO-DBO (d)CPO | `fig06_convergence_4panels.png` + per-panel `fig06_panel_a..d_*.jpg` | `fig06_panel_a_dbo.fig`, `fig06_panel_b_ipso.fig`, `fig06_panel_c_cpo_dbo.fig`, `fig06_panel_d_cpo.fig` |
| 图8  Performance boxplots, 20 instances (a-d) | `fig08_performance_boxplots_20instances.png` | `fig08_performance_boxplots.fig` |
| 图10 ECDF of optimal exposure cost | `fig10_ecdf_exposure_cost.jpg` | `fig10_ecdf.fig` |
| 图11 Worst-case risk (P90/P95/Max) | `fig11_p90_p95_max_exposure_cost.jpg` | `fig11_p90_p95_max.fig` |
| 图14 Five representative seed convergence curves | `fig14_top5_seed_convergence.jpg` | `fig14_top5_seed_convergence.fig` |
| 图20 51-node / 7-source / 5-obstacle scene | `fig20_scene_51nodes_7sources_5obstacles.jpg` | `fig20_scene_51nodes.fig` |
| 图21 Optimal path: (a) free (b) mixed obstacle model | `fig21_free_vs_obstacle_mixed_model.png`, panels `fig21_panel_a_free.jpg`, `fig21_panel_b_obstacle_mixed.jpg` | `fig21_panel_a_free.fig`, `fig21_panel_b_obstacle_mixed.fig` |
| 图23 Obstacle-retention combinations (a)1-5 (b)2-5 (c)1,3-5 (d)3-5 | `fig23_obstacle_retention_combinations.png` + per-panel `fig23_panel_*.jpg` | `fig23_panel_a_retain_12345.fig` … `fig23_panel_d_retain_345.fig` |

Figures 5, 7, 9, 12, 13, 15–19, 22 were either already correctly labelled or
are produced by the Python analysis pipeline (`python/analyze_*.py`, outputs
under `figures/revision/` and `figures/supplementary/`).

## Known pending label fixes

Tracked in `matlab/figure_label_fixes/fix_labels.m` (run once in MATLAB,
then re-export and replace the affected files):

1. 图6 panels + combined PNG: legend still reads **"Optimal Dose"** — should
   be "Optimal exposure cost". (Y-axis label is already correct.)
2. 图4 heatmap: Y-axis label reads **"Start Npde"** — typo, should be
   "Start Node".
3. 图21 panel (a) title shows **212.364 a.u.** while Table 15 reports
   **212.35 a.u.** (= 212.346). Verify the MATLAB path-sum computation and
   make figure title and table consistent (216.78 / 1070.1 for panel (b)
   already match).
