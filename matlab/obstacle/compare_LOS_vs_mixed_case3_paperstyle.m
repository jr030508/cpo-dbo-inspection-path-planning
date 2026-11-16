function compare_LOS_vs_mixed_case3_paperstyle()
% ============================================================
% compare_LOS_vs_mixed_case3_paperstyle.m
% ------------------------------------------------------------
% 第三部分分析脚本 1：
% 无障碍精细模型 vs 有障碍混合模型 总对比（论文风格）
%
% 输入：
%   1) result_case3_LOS_refined.mat
%   2) result_case3_mixed_obstacle.mat
%
% 输出目录：
%   ./compare_LOS_vs_mixed_case3_paper/
%
% 输出内容：
%   1) summary_LOS_vs_mixed.csv
%   2) delta_LOS_vs_mixed.csv
%   3) edges_removed_by_obstacle.csv
%   4) edges_added_by_obstacle.csv
%   5) compare_convergence_paper.png
%   6) compare_global_metrics_paper.png
%   7) compare_best_tours_paper.png
%   8) compare_changed_edges_paper.png
%   9) compare_LOS_vs_mixed_case3_paper.mat
% ============================================================

clc; close all;

%% ===================== 0) 文件路径 =====================
fileLOS   = fullfile(pwd, 'result_case3_LOS_refined.mat');
fileMixed = fullfile(pwd, 'result_case3_mixed_obstacle.mat');

assert(exist(fileLOS, 'file')==2, '找不到文件: %s', fileLOS);
assert(exist(fileMixed, 'file')==2, '找不到文件: %s', fileMixed);

outDir = fullfile(pwd, 'compare_LOS_vs_mixed_case3_paper');
if ~exist(outDir, 'dir')
    mkdir(outDir);
end

fprintf('\n================= 无障碍精细模型 vs 有障碍混合模型：总对比（论文风格） =================\n');
fprintf('LOS file   : %s\n', fileLOS);
fprintf('Mixed file : %s\n', fileMixed);
fprintf('Output dir : %s\n', outDir);

%% ===================== 1) 读取数据 =====================
A = load(fileLOS);
B = load(fileMixed);

reqCommon = {'gBestPath','convergence','doseMat','lenMat','costMat', ...
             'edgeType','xCoords','yCoords','sources','sourceStrength','nodeData','epsLen', ...
             'totalDose','totalLen'};
for i = 1:numel(reqCommon)
    assert(isfield(A, reqCommon{i}), '无障碍结果缺字段: %s', reqCommon{i});
    assert(isfield(B, reqCommon{i}), '有障碍结果缺字段: %s', reqCommon{i});
end

reqMixed = {'obstacles','occ','xVec','yVec','nodeRC','radGrid','gridRes','minRad'};
for i = 1:numel(reqMixed)
    assert(isfield(B, reqMixed{i}), '有障碍结果缺字段: %s', reqMixed{i});
end

%% ===================== 2) 论文风格设置 =====================
fontName = 'SimSun';
fsAxis   = 12;
fsTitle  = 14;
fsTick   = 11;
lwMain   = 2.2;
lwMinor  = 1.8;

col_free    = [0.10, 0.35, 0.80];  % 蓝
col_obs     = [0.85, 0.33, 0.10];  % 橙红
col_los     = [0.10, 0.35, 0.80];
col_astar   = [0.80, 0.15, 0.15];
col_obsFill = [0.35, 0.35, 0.35];

%% ===================== 3) 基础变量 =====================
tour_free = A.gBestPath(:).';
tour_obs  = B.gBestPath(:).';

conv_free = A.convergence(:);
conv_obs  = B.convergence(:);

doseMat_free = A.doseMat;
doseMat_obs  = B.doseMat;

lenMat_free = A.lenMat;
lenMat_obs  = B.lenMat;

costMat_free = A.costMat;
costMat_obs  = B.costMat;

edgeType_free = A.edgeType;
edgeType_obs  = B.edgeType;

xCoords = A.xCoords(:);
yCoords = A.yCoords(:);
sources = A.sources;
sourceStrength = A.sourceStrength;

obstacles = B.obstacles;
epsLen_free = A.epsLen;
epsLen_obs  = B.epsLen;

% 下面这些用于真实A*绕障绘图
occ       = B.occ;
xVec_m    = B.xVec;
yVec_m    = B.yVec;
nodeRC    = B.nodeRC;
radGrid_m = B.radGrid;
gridRes   = B.gridRes;
minRad    = B.minRad;

xMin = 0; xMax = 200;
yMin = 0; yMax = 200;

%% ===================== 4) 汇总指标 =====================
dose_free = pathSum(tour_free, doseMat_free);
dose_obs  = pathSum(tour_obs,  doseMat_obs);

len_free  = pathSum(tour_free, lenMat_free);
len_obs   = pathSum(tour_obs,  lenMat_obs);

cost_free = pathSum(tour_free, costMat_free);
cost_obs  = pathSum(tour_obs,  costMat_obs);

[firstHit_free, lastImp_free] = computeHitImprove(conv_free, conv_free(end));
[firstHit_obs,  lastImp_obs ] = computeHitImprove(conv_obs,  conv_obs(end));

nLOSEdges_free = sumEdgeTypeCount(tour_free, edgeType_free, 0);
nAStar_free    = sumEdgeTypeCount(tour_free, edgeType_free, 1);

nLOSEdges_obs  = sumEdgeTypeCount(tour_obs, edgeType_obs, 0);
nAStar_obs     = sumEdgeTypeCount(tour_obs, edgeType_obs, 1);

deltaDosePct = 100 * (dose_obs - dose_free) / max(dose_free, eps);
deltaLenPct  = 100 * (len_obs  - len_free ) / max(len_free, eps);
deltaCostPct = 100 * (cost_obs - cost_free) / max(cost_free, eps);

%% ===================== 5) 导出总表 =====================
Model = ["无障碍精细模型"; "有障碍混合模型"];
TotalDose = [dose_free; dose_obs];
TotalLength = [len_free; len_obs];
TotalCost = [cost_free; cost_obs];
FirstHitIter = [firstHit_free; firstHit_obs];
LastImproveIter = [lastImp_free; lastImp_obs];
LOSEdgeCount = [nLOSEdges_free; nLOSEdges_obs];
AStarEdgeCount = [nAStar_free; nAStar_obs];

Tsummary = table(Model, TotalDose, TotalLength, TotalCost, ...
    FirstHitIter, LastImproveIter, LOSEdgeCount, AStarEdgeCount);

writetable(Tsummary, fullfile(outDir, 'summary_LOS_vs_mixed.csv'));

Metric = ["总剂量变化率(%)"; "总距离变化率(%)"; "总代价变化率(%)"];
Value  = [deltaDosePct; deltaLenPct; deltaCostPct];
Tdelta = table(Metric, Value);
writetable(Tdelta, fullfile(outDir, 'delta_LOS_vs_mixed.csv'));

fprintf('\n===== 总体对比结果 =====\n');
disp(Tsummary);
disp(Tdelta);

%% ===================== 6) 辐射背景场 =====================
visRes = 2.0;
xVec = xMin:visRes:xMax;
yVec = yMin:visRes:yMax;
[XX, YY] = meshgrid(xVec, yVec);

radGrid = zeros(size(XX));
for s = 1:size(sources,1)
    dx = XX - sources(s,1);
    dy = YY - sources(s,2);
    radGrid = radGrid + sourceStrength(s) ./ (dx.^2 + dy.^2 + 10);
end

%% ===================== 7) 图1：收敛曲线对比 =====================
fig1 = figure('Color','w','Name','收敛曲线对比（论文风格）','Position',[100,100,980,460]);
hold on; grid on; box on;

plot(conv_free, '-', 'LineWidth', lwMain, 'Color', col_free);
plot(conv_obs,  '-', 'LineWidth', lwMain, 'Color', col_obs);

xlabel('迭代次数','FontName',fontName,'FontSize',fsAxis);
ylabel('最优巡检代价（剂量 + \epsilon·距离）','FontName',fontName,'FontSize',fsAxis);
title('无障碍精细模型与有障碍混合模型的收敛曲线对比', ...
    'FontName',fontName,'FontSize',fsTitle,'FontWeight','bold');
legend({'无障碍精细模型','有障碍混合模型'}, ...
    'Location','northeast','FontName',fontName,'FontSize',fsTick);
set(gca,'FontName',fontName,'FontSize',fsTick,'LineWidth',1.2);
hold off;

saveas(fig1, fullfile(outDir, 'compare_convergence_paper.png'));

%% ===================== 8) 图2：总体指标对比 =====================
fig2 = figure('Color','w','Name','总体指标对比（论文风格）','Position',[120,120,1180,440]);

valsDose = [dose_free, dose_obs];
valsLen  = [len_free,  len_obs];
valsCost = [cost_free, cost_obs];

subplot(1,3,1);
hold on; box on; grid on;
b = bar(valsDose, 0.58);
b.FaceColor = 'flat';
b.CData = [col_free; col_obs];
ylabel('总剂量 (\muSv)','FontName',fontName,'FontSize',fsAxis);
title('总剂量对比','FontName',fontName,'FontSize',fsTitle-1,'FontWeight','bold');
set(gca,'XTick',1:2,'XTickLabel',{'无障碍精细','有障碍混合'}, ...
    'FontName',fontName,'FontSize',fsTick,'LineWidth',1.1);
text(1.5, max(valsDose)*1.03, sprintf('\\Delta = %.2f%%', deltaDosePct), ...
    'HorizontalAlignment','center','FontName',fontName,'FontSize',11);
hold off;

subplot(1,3,2);
hold on; box on; grid on;
b = bar(valsLen, 0.58);
b.FaceColor = 'flat';
b.CData = [col_free; col_obs];
ylabel('总距离 (m)','FontName',fontName,'FontSize',fsAxis);
title('总距离对比','FontName',fontName,'FontSize',fsTitle-1,'FontWeight','bold');
set(gca,'XTick',1:2,'XTickLabel',{'无障碍精细','有障碍混合'}, ...
    'FontName',fontName,'FontSize',fsTick,'LineWidth',1.1);
text(1.5, max(valsLen)*1.03, sprintf('\\Delta = %.2f%%', deltaLenPct), ...
    'HorizontalAlignment','center','FontName',fontName,'FontSize',11);
hold off;

subplot(1,3,3);
hold on; box on; grid on;
b = bar(valsCost, 0.58);
b.FaceColor = 'flat';
b.CData = [col_free; col_obs];
ylabel('总代价','FontName',fontName,'FontSize',fsAxis);
title('总代价对比','FontName',fontName,'FontSize',fsTitle-1,'FontWeight','bold');
set(gca,'XTick',1:2,'XTickLabel',{'无障碍精细','有障碍混合'}, ...
    'FontName',fontName,'FontSize',fsTick,'LineWidth',1.1);
text(1.5, max(valsCost)*1.03, sprintf('\\Delta = %.2f%%', deltaCostPct), ...
    'HorizontalAlignment','center','FontName',fontName,'FontSize',11);
hold off;

sgtitle('无障碍精细模型与有障碍混合模型的总体指标对比', ...
    'FontName',fontName,'FontSize',fsTitle,'FontWeight','bold');

saveas(fig2, fullfile(outDir, 'compare_global_metrics_paper.png'));

%% ===================== 9) 图3：最优路径并排对比 =====================
fig3 = figure('Color','w','Name','最优路径对比（论文风格）','Position',[80,80,1450,630]);

% ---- 左：无障碍精细模型 ----
subplot(1,2,1); hold on; axis equal; box on; grid on;
xlim([xMin xMax]); ylim([yMin yMax]);

contourf(XX, YY, radGrid, 28, 'LineStyle','none');
colormap(jet);
cb = colorbar;
cb.Label.String = '辐射水平（相对值）';
cb.Label.FontName = fontName;
cb.FontName = fontName;

scatter(xCoords, yCoords, 58, 'k', 'filled', 'MarkerEdgeColor','w', 'LineWidth',0.8);
scatter(sources(:,1), sources(:,2), 125, 'r', 'pentagram', 'filled');

for i = 1:numel(xCoords)
    text(xCoords(i)+2, yCoords(i)+1, sprintf('%d', i), ...
        'Color','w','FontWeight','bold', 'BackgroundColor',[0 0 0 0.35], ...
        'HorizontalAlignment','center', 'FontName',fontName, 'FontSize',9);
end

plotTourLOS(tour_free, xCoords, yCoords, lwMinor);

xlabel('X坐标','FontName',fontName,'FontSize',fsAxis);
ylabel('Y坐标','FontName',fontName,'FontSize',fsAxis);
title(sprintf('无障碍精细模型\n总剂量 = %.3f \\muSv, 总距离 = %.1f m', dose_free, len_free), ...
    'FontName',fontName,'FontSize',fsTitle-1,'FontWeight','bold');
set(gca,'FontName',fontName,'FontSize',fsTick,'LineWidth',1.1);

% ---- 右：有障碍混合模型 ----
subplot(1,2,2); hold on; axis equal; box on; grid on;
xlim([xMin xMax]); ylim([yMin yMax]);

contourf(XX, YY, radGrid, 28, 'LineStyle','none');
colormap(jet);
cb = colorbar;
cb.Label.String = '辐射水平（相对值）';
cb.Label.FontName = fontName;
cb.FontName = fontName;

drawObstacles_CN(obstacles, col_obsFill);

scatter(xCoords, yCoords, 58, 'k', 'filled', 'MarkerEdgeColor','w', 'LineWidth',0.8);
scatter(sources(:,1), sources(:,2), 125, 'r', 'pentagram', 'filled');

for i = 1:numel(xCoords)
    text(xCoords(i)+2, yCoords(i)+1, sprintf('%d', i), ...
        'Color','w','FontWeight','bold', 'BackgroundColor',[0 0 0 0.35], ...
        'HorizontalAlignment','center', 'FontName',fontName, 'FontSize',9);
end

plotTourMixedTrueAstar(tour_obs, xCoords, yCoords, edgeType_obs, ...
    nodeRC, occ, radGrid_m, gridRes, epsLen_obs, minRad, xVec_m, yVec_m, lwMinor);

xlabel('X坐标','FontName',fontName,'FontSize',fsAxis);
ylabel('Y坐标','FontName',fontName,'FontSize',fsAxis);
title(sprintf('有障碍混合模型\n总剂量 = %.3f \\muSv, 总距离 = %.1f m, A*边 = %d', ...
    dose_obs, len_obs, nAStar_obs), ...
    'FontName',fontName,'FontSize',fsTitle-1,'FontWeight','bold');
set(gca,'FontName',fontName,'FontSize',fsTick,'LineWidth',1.1);

sgtitle('无障碍精细模型与有障碍混合模型的最优路径对比', ...
    'FontName',fontName,'FontSize',fsTitle,'FontWeight','bold');

saveas(fig3, fullfile(outDir, 'compare_best_tours_paper.png'));

%% ===================== 10) 关键变化边分析 =====================
edges_free = tourToEdges(tour_free);
edges_obs  = tourToEdges(tour_obs);

edges_free_u = sort(edges_free, 2);
edges_obs_u  = sort(edges_obs,  2);

[onlyFree, onlyObs, bothEdges] = compareEdgeSets(edges_free_u, edges_obs_u); %#ok<NASGU>

% ---- 无障碍存在但障碍后消失的边 ----
n1 = size(onlyFree,1);
Type_in_ObstacleModel = zeros(n1,1);
Dose_FreeModel = zeros(n1,1);
Dose_ObsModel  = zeros(n1,1);
Len_FreeModel  = zeros(n1,1);
Len_ObsModel   = zeros(n1,1);

for i = 1:n1
    a = onlyFree(i,1); b = onlyFree(i,2);
    Type_in_ObstacleModel(i) = edgeType_obs(a,b);
    Dose_FreeModel(i) = doseMat_free(a,b);
    Dose_ObsModel(i)  = doseMat_obs(a,b);
    Len_FreeModel(i)  = lenMat_free(a,b);
    Len_ObsModel(i)   = lenMat_obs(a,b);
end

Tremoved = table(onlyFree(:,1), onlyFree(:,2), Type_in_ObstacleModel, ...
    Dose_FreeModel, Dose_ObsModel, Len_FreeModel, Len_ObsModel, ...
    'VariableNames', {'NodeA','NodeB','Type_in_ObstacleModel','Dose_Free','Dose_Obstacle','Len_Free','Len_Obstacle'});

writetable(Tremoved, fullfile(outDir, 'edges_removed_by_obstacle.csv'));

% ---- 障碍后新增的边 ----
n2 = size(onlyObs,1);
Type_in_ObstacleModel2 = zeros(n2,1);
Dose_FreeModel2 = zeros(n2,1);
Dose_ObsModel2  = zeros(n2,1);
Len_FreeModel2  = zeros(n2,1);
Len_ObsModel2   = zeros(n2,1);

for i = 1:n2
    a = onlyObs(i,1); b = onlyObs(i,2);
    Type_in_ObstacleModel2(i) = edgeType_obs(a,b);
    Dose_FreeModel2(i) = doseMat_free(a,b);
    Dose_ObsModel2(i)  = doseMat_obs(a,b);
    Len_FreeModel2(i)  = lenMat_free(a,b);
    Len_ObsModel2(i)   = lenMat_obs(a,b);
end

Tadded = table(onlyObs(:,1), onlyObs(:,2), Type_in_ObstacleModel2, ...
    Dose_FreeModel2, Dose_ObsModel2, Len_FreeModel2, Len_ObsModel2, ...
    'VariableNames', {'NodeA','NodeB','Type_in_ObstacleModel','Dose_Free','Dose_Obstacle','Len_Free','Len_Obstacle'});

writetable(Tadded, fullfile(outDir, 'edges_added_by_obstacle.csv'));

fprintf('\n===== 关键边变化 =====\n');
fprintf('无障碍存在但有障碍消失的边数: %d\n', size(onlyFree,1));
fprintf('有障碍新增的边数: %d\n', size(onlyObs,1));

%% ===================== 11) 图4：关键变化边统计图 =====================
fig4 = figure('Color','w','Name','关键变化边统计图（论文风格）','Position',[140,140,980,430]);

subplot(1,2,1);
hold on; box on; grid on;
bar([size(onlyFree,1), size(onlyObs,1)], 0.55, 'FaceColor',[0.45 0.45 0.45]);
set(gca, 'XTick',1:2, 'XTickLabel',{'被替换边','新增边'}, ...
    'FontName',fontName,'FontSize',fsTick,'LineWidth',1.1);
ylabel('边数量','FontName',fontName,'FontSize',fsAxis);
title('路径关键变化边数量统计','FontName',fontName,'FontSize',fsTitle-1,'FontWeight','bold');
hold off;

subplot(1,2,2);
hold on; box on; grid on;
vals = [nLOSEdges_obs, nAStar_obs];
b = bar(vals, 0.55);
b.FaceColor = 'flat';
b.CData = [col_los; col_astar];
set(gca, 'XTick',1:2, 'XTickLabel',{'LOS边','A*边'}, ...
    'FontName',fontName,'FontSize',fsTick,'LineWidth',1.1);
ylabel('边数量','FontName',fontName,'FontSize',fsAxis);
title('有障碍最优路径中的边类型组成','FontName',fontName,'FontSize',fsTitle-1,'FontWeight','bold');
hold off;

sgtitle('障碍引入后的路径结构变化统计', ...
    'FontName',fontName,'FontSize',fsTitle,'FontWeight','bold');

saveas(fig4, fullfile(outDir, 'compare_changed_edges_paper.png'));

%% ===================== 12) 保存 MAT =====================
save(fullfile(outDir, 'compare_LOS_vs_mixed_case3_paper.mat'), ...
    'tour_free','tour_obs', ...
    'conv_free','conv_obs', ...
    'dose_free','dose_obs', ...
    'len_free','len_obs', ...
    'cost_free','cost_obs', ...
    'firstHit_free','firstHit_obs', ...
    'lastImp_free','lastImp_obs', ...
    'nLOSEdges_free','nLOSEdges_obs', ...
    'nAStar_free','nAStar_obs', ...
    'deltaDosePct','deltaLenPct','deltaCostPct', ...
    'Tsummary','Tdelta','Tremoved','Tadded', ...
    'onlyFree','onlyObs');

fprintf('\n✅ 总对比脚本执行完成。\n');
fprintf('输出目录：%s\n', outDir);

end


%% ============================================================
%                          Functions
%% ============================================================

function s = pathSum(path, mat)
s = 0;
for i = 1:(numel(path)-1)
    s = s + mat(path(i), path(i+1));
end
end

function [firstHit, lastImp] = computeHitImprove(conv, finalBest)
epsTol = 1e-12;
hitIdx = find(conv <= finalBest + epsTol, 1, 'first');
if isempty(hitIdx), hitIdx = numel(conv); end
firstHit = hitIdx;

d = diff(conv);
impIdx = find(d < -1e-12);
if isempty(impIdx)
    lastImp = 1;
else
    lastImp = impIdx(end) + 1;
end
end

function n = sumEdgeTypeCount(tour, edgeType, t)
n = 0;
for i = 1:(numel(tour)-1)
    if edgeType(tour(i), tour(i+1)) == t
        n = n + 1;
    end
end
end

function drawObstacles_CN(obstacles, faceColor)
for k = 1:size(obstacles,1)
    r = obstacles(k,:);
    rectangle('Position', r, 'FaceColor', [faceColor, 0.65], ...
        'EdgeColor','k', 'LineWidth',1.2);
    cx = r(1)+r(3)/2;
    cy = r(2)+r(4)/2;
    text(cx, cy, sprintf('障碍%d',k), ...
        'Color','w','FontWeight','bold', ...
        'HorizontalAlignment','center', ...
        'VerticalAlignment','middle', ...
        'FontName','SimSun','FontSize',9);
end
end

function plotTourLOS(tour, xCoords, yCoords, lw)
cols = lines(numel(tour)-1);
for k = 1:(numel(tour)-1)
    a = tour(k);
    b = tour(k+1);
    plot([xCoords(a), xCoords(b)], [yCoords(a), yCoords(b)], ...
        '-', 'LineWidth', lw, 'Color', cols(mod(k-1,size(cols,1))+1,:));
end
end

function plotTourMixedTrueAstar(tour, xCoords, yCoords, edgeType, ...
    nodeRC, occ, radGrid, gridRes, epsLen, minRad, xVec, yVec, lw)

cols = lines(numel(tour)-1);

for k = 1:(numel(tour)-1)
    a = tour(k);
    b = tour(k+1);
    c = cols(mod(k-1,size(cols,1))+1,:);

    if edgeType(a,b) == 0
        % LOS 直线边
        plot([xCoords(a), xCoords(b)], [yCoords(a), yCoords(b)], ...
            '-', 'LineWidth', lw, 'Color', c);
    else
        % A* 真实绕障边
        [pathRC, ~, ~] = astar_grid_dose_compare(nodeRC(a,:), nodeRC(b,:), ...
            occ, radGrid, gridRes, epsLen, minRad);

        pathRC2 = simplifyCollinear_compare(pathRC);
        pathRC3 = shortcutLineOfSight_compare(pathRC2, occ);

        [px, py] = gridPathToWorld_compare(pathRC3, xVec, yVec);

        plot(px, py, '--', 'LineWidth', lw+0.2, 'Color', c);
    end
end
end

function [px, py] = gridPathToWorld_compare(pathRC, xVec, yVec)
py = yVec(pathRC(:,1));
px = xVec(pathRC(:,2));
end

function [pathRC, pathLen, pathDose] = astar_grid_dose_compare(startRC, goalRC, occ, radGrid, gridRes, epsLen, minRad)

moves = [
    -1,  0, 1;
     1,  0, 1;
     0, -1, 1;
     0,  1, 1;
    -1, -1, sqrt(2);
    -1,  1, sqrt(2);
     1, -1, sqrt(2);
     1,  1, sqrt(2);
];

sr = startRC(1); sc = startRC(2);
gr = goalRC(1);  gc = goalRC(2);

gScore = inf(size(occ));
fScore = inf(size(occ));
cameR = zeros(size(occ), 'int16');
cameC = zeros(size(occ), 'int16');
inOpen = false(size(occ));
closed = false(size(occ));

gScore(sr,sc) = 0;

hPerMeter = (minRad/100) + epsLen;
fScore(sr,sc) = hPerMeter * (gridRes * hypot(gr-sr, gc-sc));
inOpen(sr,sc) = true;

while any(inOpen(:))
    [cr, cc] = argminOpen_compare(fScore, inOpen);
    if cr==gr && cc==gc
        break;
    end

    inOpen(cr,cc) = false;
    closed(cr,cc) = true;

    for k = 1:size(moves,1)
        nr = cr + moves(k,1);
        nc = cc + moves(k,2);

        if nr<1 || nr>size(occ,1) || nc<1 || nc>size(occ,2)
            continue;
        end
        if occ(nr,nc) || closed(nr,nc)
            continue;
        end

        stepLen  = moves(k,3) * gridRes;
        stepDose = stepLen * radGrid(nr,nc) / 100;
        stepCost = stepDose + epsLen * stepLen;

        tentativeG = gScore(cr,cc) + stepCost;

        if tentativeG < gScore(nr,nc)
            cameR(nr,nc) = cr;
            cameC(nr,nc) = cc;
            gScore(nr,nc) = tentativeG;

            h = hPerMeter * (gridRes * hypot(gr-nr, gc-nc));
            fScore(nr,nc) = tentativeG + h;
            inOpen(nr,nc) = true;
        end
    end
end

if ~isfinite(gScore(gr,gc))
    error('A* failed in compare script: (%d,%d)->(%d,%d) unreachable.', sr,sc,gr,gc);
end

pathRC = reconstructPath_compare(cameR, cameC, startRC, goalRC);

pathLen = 0;
pathDose = 0;
for i = 1:(size(pathRC,1)-1)
    dr = abs(pathRC(i+1,1)-pathRC(i,1));
    dc = abs(pathRC(i+1,2)-pathRC(i,2));
    if dr==1 && dc==1
        stepLen = sqrt(2)*gridRes;
    else
        stepLen = 1*gridRes;
    end
    pathLen = pathLen + stepLen;
    r = pathRC(i+1,1);
    c = pathRC(i+1,2);
    pathDose = pathDose + stepLen * radGrid(r,c) / 100;
end
end

function [r,c] = argminOpen_compare(fScore, inOpen)
tmp = fScore;
tmp(~inOpen) = inf;
[~, idx] = min(tmp(:));
[r,c] = ind2sub(size(tmp), idx);
end

function pathRC = reconstructPath_compare(cameR, cameC, startRC, goalRC)
cur = goalRC;
pathRC = cur;
while ~(cur(1)==startRC(1) && cur(2)==startRC(2))
    pr = cameR(cur(1),cur(2));
    pc = cameC(cur(1),cur(2));
    if pr==0 && pc==0
        error('Path reconstruction failed in compare script.');
    end
    cur = [double(pr), double(pc)];
    pathRC = [cur; pathRC]; %#ok<AGROW>
end
end

function path2 = simplifyCollinear_compare(pathRC)
if size(pathRC,1) <= 2
    path2 = pathRC;
    return;
end

keep = true(size(pathRC,1),1);
for i = 2:size(pathRC,1)-1
    v1 = pathRC(i,:)   - pathRC(i-1,:);
    v2 = pathRC(i+1,:) - pathRC(i,:);
    if v1(1)==v2(1) && v1(2)==v2(2)
        keep(i) = false;
    end
end
path2 = pathRC(keep,:);
end

function path3 = shortcutLineOfSight_compare(pathRC, occ)
if size(pathRC,1) <= 2
    path3 = pathRC;
    return;
end

path3 = pathRC(1,:);
i = 1;
while i < size(pathRC,1)
    j = size(pathRC,1);
    while j > i+1
        if hasLineOfSight_grid_compare(pathRC(i,:), pathRC(j,:), occ)
            break;
        end
        j = j - 1;
    end
    path3 = [path3; pathRC(j,:)]; %#ok<AGROW>
    i = j;
end
end

function ok = hasLineOfSight_grid_compare(aRC, bRC, occ)
[rr, cc] = bresenham_compare(aRC(1), aRC(2), bRC(1), bRC(2));
ok = true;
for k = 1:numel(rr)
    if occ(rr(k), cc(k))
        ok = false;
        return;
    end
end
end

function [rr, cc] = bresenham_compare(r0, c0, r1, c1)
dr = abs(r1 - r0);
dc = abs(c1 - c0);
sr = sign(r1 - r0); if sr==0, sr=1; end
sc = sign(c1 - c0); if sc==0, sc=1; end

if dc > dr
    [cc, rr] = bresenhamCore_compare(c0, r0, c1, r1, sc, sr, dc, dr);
else
    [rr, cc] = bresenhamCore_compare(r0, c0, r1, c1, sr, sc, dr, dc);
end
end

function [aa, bb] = bresenhamCore_compare(a0, b0, a1, b1, sa, sb, da, db)
aa = zeros(da+1,1);
bb = zeros(da+1,1);
err = 2*db - da;
a = a0;
b = b0;

for i = 1:(da+1)
    aa(i) = a;
    bb(i) = b;
    if err >= 0
        b = b + sb;
        err = err - 2*da;
    end
    a = a + sa;
    err = err + 2*db;
end
end

function E = tourToEdges(tour)
m = numel(tour)-1;
E = zeros(m,2);
for i = 1:m
    E(i,:) = [tour(i), tour(i+1)];
end
end

function [onlyA, onlyB, bothAB] = compareEdgeSets(A, B)
[lia, locb] = ismember(A, B, 'rows');
onlyA = A(~lia,:);
bothAB = A(lia,:);

lib = ismember(B, A, 'rows');
onlyB = B(~lib,:);

if isempty(locb) %#ok<ISMT>
end
end