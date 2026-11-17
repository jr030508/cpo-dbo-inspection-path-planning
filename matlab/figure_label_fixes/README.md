# figure_label_fixes

The authoritative, editable sources for every redrawn paper figure now live in
[`figures/fig_sources/`](../../figures/fig_sources) (MATLAB `.fig`, labels
already switched to *relative exposure cost (a.u.)*).

Three residual label issues remain; `fix_labels.m` fixes two of them
automatically and reports the third for a manual decision:

1. **图6 (convergence panels)** — legend still reads "Optimal Dose";
   the script renames it to "Optimal exposure cost" and saves patched `.fig`
   files under `output/`.
2. **图4 (heatmap)** — Y-axis label typo "Start Npde" → "Start Node".
3. **图21 panel (a)** — title shows 212.364 a.u. while the paper's Table 15
   reports 212.35 (= 212.346); the script prints the title so the
   discrepancy can be resolved deliberately (figure or table must be made
   consistent). It does **not** rewrite the number automatically.

## Usage

```matlab
cd matlab/figure_label_fixes
fix_labels
```

Patched `.fig` files appear in `output/`. Re-export each to JPG/PNG
(300 dpi, Times New Roman is already embedded in the figures) and replace the
matching files in `figures/paper/`.

> Note: an earlier iteration of this folder contained automatic
> re-plotting scripts that rebuilt 图4/6/8/10/11/14/21 from extracted or
> stored matrices. They have been retired: the author regenerated the real
> figures on the original MATLAB installation and committed the `.fig`
> sources, which supersede any re-plotting approximation.
