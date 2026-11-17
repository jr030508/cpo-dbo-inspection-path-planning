%% ============================================================
%  图14 重绘：随机初始化实验 前5名种子收敛曲线 + Partial Enlargement 放大插图
%  数据来源：从论文嵌入原图逐像素提取，平台值锚定 327.719 a.u.
%  输出：output/fig14_seed_curves_cost.png (300dpi)
%% ============================================================
scriptDir = fileparts(mfilename('fullpath'));
cd(scriptDir);
outDir = fullfile(scriptDir, 'output');
if ~exist(outDir, 'dir'); mkdir(outDir); end
FN = 'Times New Roman';

S = load(fullfile(scriptDir, 'data', 'fig14_curves.mat'));
C = S.curves5;               % 5 x 5000
leg = S.legend5;
cols = [0 0.4470 0.7410; 0.8500 0.3250 0.0980; 0.9290 0.6940 0.1250; ...
        0.4940 0.1840 0.5560; 0.4660 0.6740 0.3010];
x = (1:5000)';

f = figure('Color','w','Position',[60,60,1080,560],'Visible','off');
axMain = axes; hold(axMain, 'on'); box(axMain, 'on');
for k = 1:5
    plot(axMain, x, C(k,:), '-', 'LineWidth', 1.6, 'Color', cols(k,:));
end
xlim([0 5000]); ylim([200 1800]);
set(axMain, 'XTick', 0:500:5000, 'YTick', 200:200:1800, ...
    'FontName', FN, 'FontSize', 10, 'LineWidth', 0.9);
xlabel(axMain, 'Iterations', 'FontName', FN, 'FontSize', 12);
ylabel(axMain, 'Optimal exposure cost (a.u.)', 'FontName', FN, 'FontSize', 12);
legend(axMain, leg, 'FontName', FN, 'FontSize', 9, 'Location', 'northeast');

% 放大区域标记框（主图坐标 0-600, 300-450）
rectangle(axMain, 'Position', [0, 300, 600, 150], 'EdgeColor', 'k', 'LineWidth', 1.1);

% 插图
axIn = axes('Position', [0.24, 0.30, 0.30, 0.42]);
hold(axIn, 'on'); box(axIn, 'on');
for k = 1:5
    plot(axIn, x, C(k,:), '-', 'LineWidth', 1.5, 'Color', cols(k,:));
end
xlim(axIn, [0 600]); ylim(axIn, [300 450]);
set(axIn, 'XTick', 0:100:600, 'YTick', 300:50:450, ...
    'FontName', FN, 'FontSize', 9, 'LineWidth', 0.9);
title(axIn, 'Partial Enlargement', 'FontName', FN, 'FontSize', 11);

% 连接线（矩形角 -> 插图角，近似原图的放大指示线）
annotation('line', [0.145 0.24], [0.42 0.62], 'Color', 'k');
annotation('line', [0.145 0.24], [0.24 0.40], 'Color', 'k');

print(f, '-dpng', '-r300', fullfile(outDir, 'fig14_seed_curves_cost.png'));
close(f);
fprintf('[图14] 完成 -> output/fig14_seed_curves_cost.png\n');
