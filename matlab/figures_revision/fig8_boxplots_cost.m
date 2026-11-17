%% ============================================================
%  图8 重绘：各算法 20 扰动实例性能箱线图（4面板）
%  (a)相对通行暴露代价 (b)首次达到终值的迭代步数 (c)墙钟时间 (d)路径总距离
%  数据来源：仓库 bench results .mat（与论文同一批数据）
%  输出：output/fig8_boxplots_cost.png (300dpi)
%  注：不依赖 Statistics Toolbox，箱线图手工绘制
%% ============================================================
scriptDir = fileparts(mfilename('fullpath'));
cd(scriptDir);
outDir = fullfile(scriptDir, 'output');
if ~exist(outDir, 'dir'); mkdir(outDir); end
FN = 'Times New Roman';

algs   = {'CPO','DBO','CPO-DBO','IPSO','GA'};
files  = {'results_CPO_20instances.mat','results_DBO_20instances.mat', ...
          'results_CPO_DBO_20instances.mat','results_IPSO_20instances.mat', ...
          'results_GA_20instances.mat'};
bdir   = fullfile(scriptDir, 'data', 'bench_results');
cols   = [0 0.4470 0.7410; 0.8500 0.3250 0.0980; 0.9290 0.6940 0.1250; ...
          0.4940 0.1840 0.5560; 0.4660 0.6740 0.3010];

D = cell(1,5); L = cell(1,5); T = cell(1,5); LEN = cell(1,5);
for k = 1:5
    S = load(fullfile(bdir, files{k}));
    D{k}   = S.bestDoseRuns(:);
    L{k}   = S.firstHitIter(:);
    T{k}   = S.timeRuns(:);
    LEN{k} = S.bestLenRuns(:);
end

panelData = {D, L, T, LEN};
panelYl   = {[0 500], [0 5200], [0 500], [0 1400]};
panelYl2  = {'Optimal total exposure cost (a.u.)', 'Iterations', ...
             'Running time (s)', 'Optimal path total distance'};
panelCap  = {'(a)','(b)','(c)','(d)'};

f = figure('Color','w','Position',[60,60,1080,590],'Visible','off');
for p = 1:4
    subplot(2,2,p); hold on; box on;
    data = panelData{p};
    for k = 1:5
        drawBoxManual(data{k}, k, cols(k,:), FN);
    end
    xlim([0.4 5.6]); ylim(panelYl{p});
    set(gca, 'XTick', 1:5, 'XTickLabel', algs, 'FontName', FN, 'FontSize', 9, ...
        'LineWidth', 0.9);
    ylabel(panelYl2{p}, 'FontName', FN, 'FontSize', 10.5);
    text('Units','normalized', 'Position',[-0.12 1.04 0], 'String', panelCap{p}, ...
        'FontName', FN, 'FontSize', 13, 'FontWeight', 'bold');
end

print(f, '-dpng', '-r300', fullfile(outDir, 'fig8_boxplots_cost.png'));
close(f);
fprintf('[图8] 完成 -> output/fig8_boxplots_cost.png\n');

%% ---------- 本地函数 ----------
function drawBoxManual(x, cx, col, FN)
x = x(:);
q = manualPrctile(x, [25 50 75]);
w = 0.30;
iqr = q(3) - q(1);
lo = x(x >= q(1) - 1.5*iqr); hi = x(x <= q(3) + 1.5*iqr);
wlo = min(lo); whi = max(hi);
% 散点（带抖动）
jit = (rand(size(x))*0.16 - 0.08);
plot(cx + jit, x, 'o', 'MarkerSize', 3.2, 'MarkerFaceColor', col, ...
    'MarkerEdgeColor', col);
% 箱体
fill([cx-w cx+w cx+w cx-w], [q(1) q(1) q(3) q(3)], [0.93 0.93 0.93], ...
    'EdgeColor', [0.35 0.35 0.35], 'LineWidth', 1.0);
% 中位数
plot([cx-w cx+w], [q(2) q(2)], '-', 'Color', [0.85 0.15 0.15], 'LineWidth', 1.8);
% 须线
plot([cx cx], [q(3) whi], '-',  'Color', [0.45 0.45 0.45], 'LineWidth', 1.0);
plot([cx cx], [q(1) wlo], '-',  'Color', [0.45 0.45 0.45], 'LineWidth', 1.0);
plot([cx-0.14 cx+0.14], [whi whi], '-', 'Color', [0.45 0.45 0.45], 'LineWidth', 1.0);
plot([cx-0.14 cx+0.14], [wlo wlo], '-', 'Color', [0.45 0.45 0.45], 'LineWidth', 1.0);
end

function q = manualPrctile(x, pcts)
% MATLAB prctile 兼容实现（线性插值）
x = sort(x(:));
n = numel(x);
q = zeros(size(pcts));
for i = 1:numel(pcts)
    pos = pcts(i)/100 * (n - 1) + 1;
    lo = floor(pos); hi = ceil(pos);
    if lo < 1; lo = 1; end
    if hi > n; hi = n; end
    t = pos - lo;
    q(i) = x(lo)*(1-t) + x(hi)*t;
end
end
