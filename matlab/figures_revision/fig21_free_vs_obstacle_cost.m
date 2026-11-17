%% ============================================================
%  图21 重绘：无障碍精细模型 vs 有障碍混合模型 最优路径对比
%  (a)无障碍精细模型  (b)有障碍混合模型（含A*绕障边虚线）
%  数据来源：原始 result_case3_*.mat（pathSum 校验 = 表15 数值）
%  输出：output/fig21_free_vs_obstacle_cost.png (300dpi)
%% ============================================================
scriptDir = fileparts(mfilename('fullpath'));
cd(scriptDir);
outDir = fullfile(scriptDir, 'output');
if ~exist(outDir, 'dir'); mkdir(outDir); end
FN = 'Times New Roman';

A = load(fullfile(scriptDir, 'data', 'fig21_LOS.mat'));    % 无障碍精细模型
B = load(fullfile(scriptDir, 'data', 'fig21_MIXED.mat'));  % 有障碍混合模型

tour_free = A.gBestPath(:).';
tour_obs  = B.gBestPath(:).';
xCoords   = A.xCoords(:);
yCoords   = A.yCoords(:);
sources   = A.sources;
sourceStrength = A.sourceStrength;
obstacles = B.obstacles;

% 与表15 一致的汇总数值（pathSum 复核）
dose_free = pathSum(tour_free, A.doseMat);
len_free  = pathSum(tour_free, A.lenMat);
dose_obs  = pathSum(tour_obs,  B.doseMat);
len_obs   = pathSum(tour_obs,  B.lenMat);
nAStar    = sumEdgeType(tour_obs, B.edgeType, 1);
fprintf('[图21] 复核: 无障碍 %.2f a.u. / %.1f m; 有障碍 %.2f a.u. / %.1f m; A*边 %d\n', ...
    dose_free, len_free, dose_obs, len_obs, nAStar);

%% 辐射背景场
visRes = 2.0; xMin = 0; xMax = 200; yMin = 0; yMax = 200;
xVec = xMin:visRes:xMax;  yVec = yMin:visRes:yMax;
[XX, YY] = meshgrid(xVec, yVec);
radGrid = zeros(size(XX));
for s = 1:size(sources,1)
    dx = XX - sources(s,1);  dy = YY - sources(s,2);
    radGrid = radGrid + sourceStrength(s) ./ (dx.^2 + dy.^2 + 10);
end

lwMinor = 1.6;

f = figure('Color','w','Position',[60,60,1180,640],'Visible','off');

% ---------------- (a) 无障碍精细模型 ----------------
subplot(1,2,1); hold on; axis equal; box on; grid on;
set(gca, 'Position', [0.06, 0.08, 0.40, 0.76]);
text('Units','normalized','Position',[-0.02,1.05,0],'String','(a)', ...
    'FontName',FN,'FontSize',13,'FontWeight','bold');
xlim([xMin xMax]); ylim([yMin yMax]);
contourf(XX, YY, radGrid, 28, 'LineStyle', 'none');
colormap(jet);
cb = colorbar;
cb.Label.String = 'Radiation level (relative value)';
cb.Label.FontName = FN; cb.Label.FontSize = 11;
cb.FontName = FN; cb.FontSize = 10;

scatter(xCoords, yCoords, 42, 'k', 'filled', 'MarkerEdgeColor', 'w', 'LineWidth', 0.7);
for s = 1:size(sources,1)
    scatter(sources(s,1), sources(s,2), 130, [1 0.85 0.1], 'pentagram', 'filled', ...
        'MarkerEdgeColor', [0.8 0.2 0]);
    text(sources(s,1), sources(s,2)-9, sprintf('Source %d', s), 'Color', 'r', ...
        'FontName', FN, 'FontSize', 8.5, 'HorizontalAlignment', 'center', ...
        'FontWeight', 'bold');
end
for i = 1:numel(xCoords)
    text(xCoords(i)+2, yCoords(i)+1, sprintf('%d', i), 'Color', 'w', ...
        'FontWeight','bold', 'BackgroundColor', [0 0 0 0.35], ...
        'HorizontalAlignment','center', 'FontName', FN, 'FontSize', 8);
end
plotTourLOS(tour_free, xCoords, yCoords, lwMinor);

xlabel('X', 'FontName', FN, 'FontSize', 12);
ylabel('Y', 'FontName', FN, 'FontSize', 12);
title(sprintf('Total exposure cost = %.2f a.u.; Total Distance = %.1f', ...
    dose_free, len_free), 'FontName', FN, 'FontSize', 11, 'FontWeight', 'bold');
set(gca, 'FontName', FN, 'FontSize', 9, 'LineWidth', 1.0);

% ---------------- (b) 有障碍混合模型 ----------------
subplot(1,2,2); hold on; axis equal; box on; grid on;
set(gca, 'Position', [0.55, 0.08, 0.40, 0.76]);
text('Units','normalized','Position',[-0.02,1.05,0],'String','(b)', ...
    'FontName',FN,'FontSize',13,'FontWeight','bold');
xlim([xMin xMax]); ylim([yMin yMax]);
contourf(XX, YY, radGrid, 28, 'LineStyle', 'none');
colormap(jet);
cb = colorbar;
cb.Label.String = 'Radiation level (relative value)';
cb.Label.FontName = FN; cb.Label.FontSize = 11;
cb.FontName = FN; cb.FontSize = 10;

for k = 1:size(obstacles,1)
    rectangle('Position', obstacles(k,:), 'FaceColor', [0.55 0.55 0.55 0.9], ...
        'EdgeColor', 'k', 'LineWidth', 1.1);
end
scatter(xCoords, yCoords, 42, 'k', 'filled', 'MarkerEdgeColor', 'w', 'LineWidth', 0.7);
for s = 1:size(sources,1)
    scatter(sources(s,1), sources(s,2), 130, [1 0.85 0.1], 'pentagram', 'filled', ...
        'MarkerEdgeColor', [0.8 0.2 0]);
    text(sources(s,1), sources(s,2)-9, sprintf('Source %d', s), 'Color', 'r', ...
        'FontName', FN, 'FontSize', 8.5, 'HorizontalAlignment', 'center', ...
        'FontWeight', 'bold');
end
for i = 1:numel(xCoords)
    text(xCoords(i)+2, yCoords(i)+1, sprintf('%d', i), 'Color', 'w', ...
        'FontWeight','bold', 'BackgroundColor', [0 0 0 0.35], ...
        'HorizontalAlignment','center', 'FontName', FN, 'FontSize', 8);
end
plotTourMixed(tour_obs, xCoords, yCoords, B.edgeType, B.nodeRC, B.occ, ...
    B.radGrid, B.gridRes, B.epsLen, B.minRad, B.xVec, B.yVec, lwMinor);

xlabel('X', 'FontName', FN, 'FontSize', 12);
ylabel('Y', 'FontName', FN, 'FontSize', 12);
title(sprintf('Total exposure cost = %.2f a.u.; Total Distance = %.1f', ...
    dose_obs, len_obs), 'FontName', FN, 'FontSize', 11, 'FontWeight', 'bold');
set(gca, 'FontName', FN, 'FontSize', 9, 'LineWidth', 1.0);

% 总图例（放顶部，独立隐藏坐标轴）
ah = axes('Position', [0.03, 0.905, 0.94, 0.08], 'Visible', 'off');
hold(ah, 'on');
h1 = patch(ah, [0 1 1], [0 0 1], [0.55 0.55 0.55], 'EdgeColor', 'k', 'LineWidth', 1.1);
h2 = plot(ah, NaN, NaN, 'pentagram', 'MarkerFaceColor', [1 0.85 0.1], ...
    'MarkerEdgeColor', [0.8 0.2 0], 'MarkerSize', 10, 'LineStyle', 'none');
h3 = plot(ah, NaN, NaN, '-',  'Color', [0.93 0.47 0.11], 'LineWidth', 2.2);
h4 = plot(ah, NaN, NaN, '--', 'Color', [0.20 0.70 0.30], 'LineWidth', 2.2);
h5 = plot(ah, NaN, NaN, 'o', 'MarkerFaceColor', 'k', 'MarkerEdgeColor', 'w', ...
    'MarkerSize', 6, 'LineStyle', 'none');
lg = legend(ah, [h1 h2 h3 h4 h5], {'Obstacle', 'Radiation source', 'LOS edge', ...
    'A* detour edge', 'Inspection node'}, 'FontName', FN, 'FontSize', 10.5, ...
    'Orientation', 'landscape', 'NumColumns', 5, 'Box', 'off');
set(lg, 'Units', 'normalized', 'Position', [0.02, 0.02, 0.96, 0.9]);

print(f, '-dpng', '-r300', fullfile(outDir, 'fig21_free_vs_obstacle_cost.png'));
close(f);
fprintf('[图21] 完成 -> output/fig21_free_vs_obstacle_cost.png\n');

%% ---------- 本地函数 ----------
function total = pathSum(tour, M)
t = tour(:)';
total = 0;
for i = 1:(numel(t)-1)
    total = total + M(t(i), t(i+1));
end
end

function n = sumEdgeType(tour, edgeType, t)
n = 0;
for i = 1:(numel(tour)-1)
    if edgeType(tour(i), tour(i+1)) == t
        n = n + 1;
    end
end
end

function plotTourLOS(tour, xCoords, yCoords, lw)
cols = lines(numel(tour)-1);
for k = 1:(numel(tour)-1)
    a = tour(k); b = tour(k+1);
    plot([xCoords(a), xCoords(b)], [yCoords(a), yCoords(b)], ...
        '-', 'LineWidth', lw, 'Color', cols(mod(k-1,size(cols,1))+1,:));
end
end

function plotTourMixed(tour, xCoords, yCoords, edgeType, nodeRC, occ, ...
    radGrid, gridRes, epsLen, minRad, xVec, yVec, lw)
cols = lines(numel(tour)-1);
for k = 1:(numel(tour)-1)
    a = tour(k); b = tour(k+1);
    c = cols(mod(k-1,size(cols,1))+1,:);
    if edgeType(a,b) == 0
        plot([xCoords(a), xCoords(b)], [yCoords(a), yCoords(b)], ...
            '-', 'LineWidth', lw, 'Color', c);
    else
        [pathRC, ~] = astar_grid(nodeRC(a,:), nodeRC(b,:), occ, radGrid, gridRes, epsLen, minRad);
        pathRC2 = simplifyCollinear(pathRC);
        pathRC3 = shortcutLOS(pathRC2, occ);
        [px, py] = toWorld(pathRC3, xVec, yVec);
        plot(px, py, '--', 'LineWidth', lw+0.2, 'Color', c);
    end
end
end

function [px, py] = toWorld(pathRC, xVec, yVec)
py = yVec(pathRC(:,1));
px = xVec(pathRC(:,2));
end

function [pathRC, pathLen] = astar_grid(startRC, goalRC, occ, radGrid, gridRes, epsLen, minRad)
moves = [ -1 0 1; 1 0 1; 0 -1 1; 0 1 1; -1 -1 sqrt(2); -1 1 sqrt(2); 1 -1 sqrt(2); 1 1 sqrt(2) ];
sr = startRC(1); sc = startRC(2); gr = goalRC(1); gc = goalRC(2);
gScore = inf(size(occ)); fScore = inf(size(occ));
cameR = zeros(size(occ),'int16'); cameC = zeros(size(occ),'int16');
inOpen = false(size(occ)); closed = false(size(occ));
gScore(sr,sc) = 0;
hPerMeter = (minRad/100) + epsLen;
fScore(sr,sc) = hPerMeter * (gridRes * hypot(gr-sr, gc-sc));
inOpen(sr,sc) = true;
while any(inOpen(:))
    tmp = fScore; tmp(~inOpen) = inf;
    [~, idx] = min(tmp(:));
    [cr, cc] = ind2sub(size(tmp), idx);
    if cr==gr && cc==gc, break; end
    inOpen(cr,cc) = false; closed(cr,cc) = true;
    for k = 1:size(moves,1)
        nr = cr + moves(k,1); nc = cc + moves(k,2);
        if nr<1 || nr>size(occ,1) || nc<1 || nc>size(occ,2), continue; end
        if occ(nr,nc) || closed(nr,nc), continue; end
        stepLen  = moves(k,3) * gridRes;
        stepDose = stepLen * radGrid(nr,nc) / 100;
        stepCost = stepDose + epsLen * stepLen;
        tentativeG = gScore(cr,cc) + stepCost;
        if tentativeG < gScore(nr,nc)
            cameR(nr,nc) = cr; cameC(nr,nc) = cc;
            gScore(nr,nc) = tentativeG;
            fScore(nr,nc) = tentativeG + hPerMeter * (gridRes * hypot(gr-nr, gc-nc));
            inOpen(nr,nc) = true;
        end
    end
end
% 回溯
cur = [gr, gc]; pathRC = cur;
while ~(cur(1)==sr && cur(2)==sc)
    pr = cameR(cur(1),cur(2)); pc = cameC(cur(1),cur(2));
    if pr==0 && pc==0, break; end
    cur = [double(pr), double(pc)];
    pathRC = [cur; pathRC]; %#ok<AGROW>
end
% 计算长度
pathLen = 0;
for i = 1:size(pathRC,1)-1
    stepLen = hypot(pathRC(i+1,1)-pathRC(i,1), pathRC(i+1,2)-pathRC(i,2)) * gridRes;
    pathLen = pathLen + stepLen;
end
end

function path2 = simplifyCollinear(pathRC)
if size(pathRC,1) <= 2, path2 = pathRC; return; end
keep = true(size(pathRC,1),1);
for i = 2:size(pathRC,1)-1
    v1 = pathRC(i,:) - pathRC(i-1,:);
    v2 = pathRC(i+1,:) - pathRC(i,:);
    if v1(1)==v2(1) && v1(2)==v2(2), keep(i) = false; end
end
path2 = pathRC(keep,:);
end

function path3 = shortcutLOS(pathRC, occ)
if size(pathRC,1) <= 2, path3 = pathRC; return; end
path3 = pathRC(1,:);
i = 1;
while i < size(pathRC,1)
    j = size(pathRC,1);
    while j > i+1
        if hasLOS(pathRC(i,:), pathRC(j,:), occ), break; end
        j = j - 1;
    end
    path3 = [path3; pathRC(j,:)]; %#ok<AGROW>
    i = j;
end
end

function ok = hasLOS(aRC, bRC, occ)
[rr, cc] = bres(aRC(1), aRC(2), bRC(1), bRC(2));
ok = true;
for k = 1:numel(rr)
    if occ(rr(k), cc(k)), ok = false; return; end
end
end

function [rr, cc] = bres(r0, c0, r1, c1)
dr = abs(r1-r0); dc = abs(c1-c0);
sr = sign(r1-r0); if sr==0, sr=1; end
sc = sign(c1-c0); if sc==0, sc=1; end
if dc > dr
    [cc, rr] = bresCore(c0, r0, c1, r1, sc, sr, dc, dr);
else
    [rr, cc] = bresCore(r0, c0, r1, c1, sr, sc, dr, dc);
end
end

function [aa, bb] = bresCore(a0, b0, a1, b1, sa, sb, da, db)
aa = zeros(da+1,1); bb = zeros(da+1,1);
err = 2*db - da; a = a0; b = b0;
for i = 1:(da+1)
    aa(i) = a; bb(i) = b;
    if err >= 0
        b = b + sb; err = err - 2*da;
    end
    a = a + sa; err = err + 2*db;
end
end
