%% ============================================================
%  图4 重绘：30节点对间暴露代价矩阵热力图
%  原图为 "有效剂量(μSv)" 标注，本版统一为相对暴露代价 (a.u.)
%  数据来源：原始 30热力图.fig（与论文完全一致），备用 doseMatrix30
%  输出：output/fig4_30node_cost_heatmap.png (300dpi)
%% ============================================================
scriptDir = fileparts(mfilename('fullpath'));
cd(scriptDir);
outDir = fullfile(scriptDir, 'output');
if ~exist(outDir, 'dir'); mkdir(outDir); end
FN = 'Times New Roman';

%% ---- 1) 读取原始矩阵（优先从 .fig 提取，保证与论文一致）----
M = [];
figPath = fullfile(scriptDir, 'data', 'fig4_30heatmap_original.fig');
if exist(figPath, 'file')
    try
        f0 = openfig(figPath, 'invisible');
        ax0 = get(f0, 'CurrentAxes');
        im = findobj(ax0, 'Type', 'image');
        if isempty(im)
            im = findobj(ax0, 'Type', 'surface');
        end
        M = get(im(1), 'CData');
        close(f0);
        fprintf('[图4] 已从原始 .fig 提取矩阵: %d x %d\n', size(M,1), size(M,2));
    catch ME
        fprintf('[图4] .fig 读取失败(%s)，改用备用矩阵\n', ME.message);
        M = [];
    end
end
if isempty(M)
    S = load(fullfile(scriptDir, 'data', 'fig4_dose30.mat'));
    M = S.doseMatrix30;
end
n = size(M, 1);
M = double(M);

%% ---- 2) 热点框（前10大非对角单元对，运行时计算）----
tmp = M; tmp(1:n+1:end) = -inf;
vals = [];
pr = []; pc = [];
for i = 1:n
    for j = i+1:n
        vals(end+1) = tmp(i,j); %#ok<*AGROW>
        pr(end+1) = i; pc(end+1) = j;
    end
end
[~, order] = sort(vals, 'descend');
hotList = order(1:10);

%% ---- 3) 绘图 ----
f = figure('Color','w','Position',[60,60,1150,620],'Visible','off');
imagesc(M);
colormap(jet);
cb = colorbar;
cb.Label.String = 'Effective exposure cost (a.u.)';
cb.Label.FontName = FN;
cb.Label.FontSize = 12;
cb.FontName = FN;
cb.FontSize = 10;

set(gca, 'XTick', 1:n, 'YTick', 1:n, 'FontName', FN, 'FontSize', 8, ...
    'DataAspectRatio', [1 1 1], 'LineWidth', 0.8);
xlabel('Target Node', 'FontName', FN, 'FontSize', 13);
ylabel('Start Node',  'FontName', FN, 'FontSize', 13);
hold on;

% 热点红虚线框（对称位置都画）
for k = 1:numel(hotList)
    idx = hotList(k);
    i = pr(idx); j = pc(idx);
    rectangle('Position', [j-0.5, i-0.5, 1, 1], 'EdgeColor','r', ...
        'LineStyle','--', 'LineWidth', 1.3);
    rectangle('Position', [i-0.5, j-0.5, 1, 1], 'EdgeColor','r', ...
        'LineStyle','--', 'LineWidth', 1.3);
end

% 数值标注（与原图一致的阈值）
labelThr = 16.0;
for i = 1:n
    for j = 1:n
        if i ~= j && M(i,j) >= labelThr
            text(j, i, sprintf('%.0f', M(i,j)), 'Color', [0.55 0 0], ...
                'FontSize', 7, 'FontName', FN, ...
                'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle');
        end
    end
end

print(f, '-dpng', '-r300', fullfile(outDir, 'fig4_30node_cost_heatmap.png'));
close(f);
fprintf('[图4] 完成 -> output/fig4_30node_cost_heatmap.png\n');
