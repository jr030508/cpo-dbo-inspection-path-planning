% FIX_LABELS  One-shot label fixes for the redrawn paper figures.
% Run this ONCE on the machine that holds the original .fig workflow, then
% re-export the affected JPG/PNG files (File > Export, or print -djpeg/-dpng)
% and replace the matching files in figures/paper/.
%
% Fixes (see figures/README.md):
%   1) fig06 panels      : legend "Optimal Dose"      -> "Optimal exposure cost"
%   2) fig04 heatmap     : Y label "Start Npde"       -> "Start Node"
%   3) fig21 panel (a)   : title value 212.364 a.u.   -> verify against Table 15
%                          (paper reports 212.35 = 212.346); fix the source
%                          number in the title string or the table, they must
%                          agree.

HERE = fileparts(mfilename('fullpath'));
SRC = fullfile(HERE, '..', '..', 'figures', 'fig_sources');
OUT = fullfile(HERE, 'output');
if ~exist(OUT, 'dir'); mkdir(OUT); end

%% --- 1) fig06: rename legend text in every convergence panel -------------
panels = {'fig06_panel_a_dbo.fig', 'fig06_panel_b_ipso.fig', ...
          'fig06_panel_c_cpo_dbo.fig', 'fig06_panel_d_cpo.fig'};
for k = 1:numel(panels)
    f = openfig(fullfile(SRC, panels{k}), 'visible', 'off');
    lines = findobj(f, 'Type', 'line');
    for h = lines'
        if isequal(get(h, 'DisplayName'), 'Optimal Dose')
            set(h, 'DisplayName', 'Optimal exposure cost');
        end
    end
    lg = findobj(f, 'Type', 'legend');
    if ~isempty(lg); legend(lg(1), 'refresh'); end
    savefig(f, fullfile(OUT, panels{k}));           % patched .fig
    close(f);
    fprintf('fig06 legend fixed : %s\n', panels{k});
end

%% --- 2) fig04: fix Y-axis typo -------------------------------------------
f = openfig(fullfile(SRC, 'fig04_exposure_cost_heatmap.fig'), 'visible', 'off');
ax = findobj(f, 'Type', 'axes');
for a = ax'
    yl = get(a, 'YLabel');
    if isequal(get(yl, 'String'), 'Start Npde')
        set(yl, 'String', 'Start Node');
    end
end
savefig(f, fullfile(OUT, 'fig04_exposure_cost_heatmap.fig'));
close(f);
fprintf('fig04 Y label fixed : Start Npde -> Start Node\n');

%% --- 3) fig21 panel (a): reported value check -----------------------------
% The title currently reads 212.364 a.u.; Table 15 reports 212.35 (=212.346).
% Open the source and print the title so the discrepancy can be resolved by
% hand (either recompute pathSum with the paper's edge-cost formula, or fix
% the table). Automatic rewrite is intentionally NOT done.
f = openfig(fullfile(SRC, 'fig21_panel_a_free.fig'), 'visible', 'off');
ttl = get(get(f, 'CurrentAxes'), 'Title');
fprintf('fig21(a) current title: %s\n', get(ttl, 'String'));
fprintf('  -> paper Table 15 expects 212.35 / 1046.5. Resolve manually.\n');
close(f);

fprintf('\nDone. Patched .fig files are in:\n  %s\n', OUT);
fprintf('Re-export JPG/PNG at 300 dpi (print -djpeg -r300 / -dpng -r300)\n');
fprintf('and replace the corresponding files in figures/paper/.\n');
