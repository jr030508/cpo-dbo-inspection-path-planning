%% ============================================================
%  图10 重绘：各算法最优总暴露代价 ECDF 曲线（20个扰动实例）
%  输出：output/fig10_ecdf_cost.png (300dpi)
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
cols  = [0 0.4470 0.7410; 0.8500 0.3250 0.0980; 0.9290 0.6940 0.1250; ...
         0.4940 0.1840 0.5560; 0.4660 0.6740 0.3010];

f = figure('Color','w','Position',[60,60,1080,560],'Visible','off');
hold on; box on;
for k = 1:5
    S = load(fullfile(bdir, files{k}));
    x = sort(S.bestDoseRuns(:));
    n = numel(x);
    F = ((1:n)' - 0.5) / n;
    X = [x(1)-8; x; x(end)+8];
    Y = [0; F; 1];
    plot(X, Y, '-', 'LineWidth', 1.7, 'Color', cols(k,:));
end
xlim([320 460]); ylim([0 1]);
set(gca, 'XTick', 320:20:460, 'YTick', 0:0.1:1, 'FontName', FN, 'FontSize', 10, ...
    'LineWidth', 0.9);
xlabel('Optimal exposure cost (a.u.)', 'FontName', FN, 'FontSize', 12);
ylabel('ECDF', 'FontName', FN, 'FontSize', 12);
legend(algs, 'FontName', FN, 'FontSize', 9, 'Location', 'southeast');

print(f, '-dpng', '-r300', fullfile(outDir, 'fig10_ecdf_cost.png'));
close(f);
fprintf('[图10] 完成 -> output/fig10_ecdf_cost.png\n');
