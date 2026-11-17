%% ============================================================
%  图11 重绘：最差情形风险对比 P90 / P95 / Max（20个扰动实例）
%  输出：output/fig11_risk_bars_cost.png (300dpi)
%% ============================================================
scriptDir = fileparts(mfilename('fullpath'));
cd(scriptDir);
outDir = fullfile(scriptDir, 'output');
if ~exist(outDir, 'dir'); mkdir(outDir); end
FN = 'Times New Roman';

algs  = {'CPO','DBO','CPO-DBO','IPSO','GA'};
files = {'results_CPO_20instances.mat','results_DBO_20instances.mat', ...
         'results_CPO_DBO_20instances.mat','results_IPSO_20instances.mat', ...
         'results_GA_20instances.mat'};
bdir  = fullfile(scriptDir, 'data', 'bench_results');

P90 = zeros(1,5); P95 = zeros(1,5); MX = zeros(1,5);
for k = 1:5
    S = load(fullfile(bdir, files{k}));
    x = sort(S.bestDoseRuns(:));
    P90(k) = manualPrctile(x, 90);
    P95(k) = manualPrctile(x, 95);
    MX(k)  = max(x);
end

D = [P90; P95; MX];
cols = [0 0.4470 0.7410; 0.8500 0.3250 0.0980; 0.9290 0.6940 0.1250];

f = figure('Color','w','Position',[60,60,1080,550],'Visible','off');
b = bar(D, 'grouped', 'Width', 0.85);
for k = 1:3
    b(k).FaceColor = cols(k,:);
    b(k).EdgeColor = 'none';
end
xlim([0.35 5.65]); ylim([0 500]);
set(gca, 'XTick', 1:5, 'XTickLabel', algs, 'FontName', FN, 'FontSize', 10, ...
    'LineWidth', 0.9);
ylabel('Optimal exposure cost (a.u.)', 'FontName', FN, 'FontSize', 12);
legend({'P90','P95','Max'}, 'FontName', FN, 'FontSize', 10, 'Location', 'northwest');

print(f, '-dpng', '-r300', fullfile(outDir, 'fig11_risk_bars_cost.png'));
close(f);
fprintf('[图11] 完成 -> output/fig11_risk_bars_cost.png\n');

function q = manualPrctile(x, p)
x = sort(x(:));
n = numel(x);
pos = p/100 * (n - 1) + 1;
lo = max(1, floor(pos)); hi = min(n, ceil(pos));
t = pos - lo;
q = x(lo)*(1-t) + x(hi)*t;
end
