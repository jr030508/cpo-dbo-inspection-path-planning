function analyze_dose_breakdown_case3_paperstyle_v2()
% ============================================================
% analyze_dose_breakdown_case3_paperstyle_v2.m
% ------------------------------------------------------------
% 第三部分分析脚本 3（升级版，修正端点对齐）：
% 剂量分解分析（论文风格）
%
% 改进点：
%   1) 关键A*边局部机制图改为“单独成图”
%   2) 图中明确标注障碍编号和节点编号
%   3) A*绕障路径首尾点强制对齐到真实节点坐标
%
% 目标：
%   1) 分析有障碍混合模型中哪些边是高剂量贡献边
%   2) 分析关键A*边相较无障碍直线边，是降剂量还是主要避障
%   3) 分析障碍引入后发生替换的关键边
%
% 输入：
%   1) result_case3_LOS_refined.mat
%   2) result_case3_mixed_obstacle.mat
%
% 输出目录：
%   ./analyze_dose_breakdown_case3_paper_v2/
% ============================================================

clc; close all;

%% ===================== 0) 文件路径 =====================
fileLOS   = fullfile(pwd, 'result_case3_LOS_refined.mat');
fileMixed = fullfile(pwd, 'result_case3_mixed_obstacle.mat');

assert(exist(fileLOS, 'file')==2, '找不到文件: %s', fileLOS);
assert(exist(fileMixed, 'file')==2, '找不到文件: %s', fileMixed);

outDir = fullfile(pwd, 'analyze_dose_breakdown_case3_paper_v2');
if ~exist(outDir, 'dir')
    mkdir(outDir);
end

fprintf('\n================= 剂量分解分析（论文风格，升级版） =================\n');
fprintf('LOS file   : %s\n', fileLOS);
fprintf('Mixed file : %s\n', fileMixed);
fprintf('Output dir : %s\n', outDir);

%% ===================== 1) 读取数据 =====================
A = load(fileLOS);
B = load(fileMixed);

reqLOS = {'gBestPath','doseMat','lenMat','costMat','edgeType','totalDose','totalLen', ...
          'xCoords','yCoords','sources','sourceStrength'};
reqMix = {'gBestPath','doseMat','lenMat','costMat','edgeType','totalDose','totalLen', ...
          'xCoords','yCoords','sources','sourceStrength', ...
          'obstacles','occ','xVec','yVec','nodeRC','radGrid','gridRes','minRad','epsLen'};

for i = 1:numel(reqLOS)
    assert(isfield(A, reqLOS{i}), '无障碍结果缺字段: %s', reqLOS{i});
end
for i = 1:numel(reqMix)
    assert(isfield(B, reqMix{i}), '有障碍结果缺字段: %s', reqMix{i});
end

%% ===================== 2) 风格参数 =====================
fontName = 'SimSun';
fsAxis   = 12;
fsTitle  = 14;
fsTick   = 11;

col_los   = [0.10, 0.35, 0.80];
col_astar = [0.80, 0.15, 0.15];
col_free  = [0.15, 0.55, 0.55];
col_gray  = [0.45, 0.45, 0.45];
col_obsFill = [0.35, 0.35, 0.35];

topK_all    = 15;   % 总体高贡献边显示前K条
topK_astar  = 8;    % 关键A*边分析前K条
topK_change = 8;    % 变化边分析前K条
localPad    = 40;   % 关键边局部图显示范围半宽度

%% ===================== 3) 提取变量 =====================
tour_free = A.gBestPath(:).';
tour_obs  = B.gBestPath(:).';

doseMat_free = A.doseMat;
lenMat_free  = A.lenMat;
edgeType_free = A.edgeType;

doseMat_obs = B.doseMat;
lenMat_obs  = B.lenMat;
edgeType_obs = B.edgeType;

totalDose_free = A.totalDose;
totalLen_free  = A.totalLen;

totalDose_obs = B.totalDose;
totalLen_obs  = B.totalLen;

xCoords = A.xCoords(:);
yCoords = A.yCoords(:);
sources = A.sources;
sourceStrength = A.sourceStrength;

obstacles = B.obstacles;
occ       = B.occ;
xVec      = B.xVec;
yVec      = B.yVec;
nodeRC    = B.nodeRC;
radGrid   = B.radGrid;
gridRes   = B.gridRes;
minRad    = B.minRad;
epsLen    = B.epsLen;

xMin = 0; xMax = 200;
yMin = 0; yMax = 200;

%% ===================== 4) 有障碍最优路径逐边分解 =====================
[doseList_obs, lenList_obs, typeList_obs, edgePairs_obs] = ...
    getPathEdgeDetail(tour_obs, doseMat_obs, lenMat_obs, edgeType_obs);

doseRatio_obs = 100 * doseList_obs / max(totalDose_obs, eps);
lenRatio_obs  = 100 * lenList_obs  / max(totalLen_obs,  eps);

EdgeIndex = (1:numel(doseList_obs)).';
NodeA = edgePairs_obs(:,1);
NodeB = edgePairs_obs(:,2);
EdgeType = strings(numel(typeList_obs),1);
EdgeType(typeList_obs==0) = "LOS";
EdgeType(typeList_obs==1) = "A*";

Tall = table(EdgeIndex, NodeA, NodeB, EdgeType, doseList_obs, lenList_obs, ...
    doseRatio_obs, lenRatio_obs, ...
    'VariableNames', {'EdgeIndex','NodeA','NodeB','EdgeType','Dose','Length','DoseRatioPct','LengthRatioPct'});

writetable(Tall, fullfile(outDir, 'dose_breakdown_mixed_all_edges.csv'));

%% ===================== 5) 高贡献边排序 =====================
[sortedDose_all, idxSort_all] = sort(doseList_obs, 'descend');
kShow_all = min(topK_all, numel(sortedDose_all));
idxTop_all = idxSort_all(1:kShow_all);

doseTop_all = doseList_obs(idxTop_all);
lenTop_all  = lenList_obs(idxTop_all);
typeTop_all = typeList_obs(idxTop_all);
pairTop_all = edgePairs_obs(idxTop_all,:);

labelsTop_all = strings(kShow_all,1);
for i = 1:kShow_all
    labelsTop_all(i) = sprintf('(%d,%d)', pairTop_all(i,1), pairTop_all(i,2));
end

TtopAll = Tall(idxTop_all,:);
writetable(TtopAll, fullfile(outDir, 'top_dose_edges_mixed.csv'));

%% ===================== 6) 关键 A* 边分析 =====================
idxAstar = find(typeList_obs == 1);

if isempty(idxAstar)
    warning('有障碍最优路径中没有 A* 边，关键A*边分析将跳过。');
    Tastar = table();
else
    [~, orderAstar] = sort(doseList_obs(idxAstar), 'descend');
    idxAstar_sorted = idxAstar(orderAstar);

    kShow_astar = min(topK_astar, numel(idxAstar_sorted));
    idxTop_astar = idxAstar_sorted(1:kShow_astar);

    EdgeIndex_astar = idxTop_astar(:);
    NodeA_astar = edgePairs_obs(idxTop_astar,1);
    NodeB_astar = edgePairs_obs(idxTop_astar,2);

    Dose_Free = zeros(kShow_astar,1);
    Dose_Obs  = zeros(kShow_astar,1);
    Len_Free  = zeros(kShow_astar,1);
    Len_Obs   = zeros(kShow_astar,1);

    DoseChangePct = zeros(kShow_astar,1);
    LenChangePct  = zeros(kShow_astar,1);

    Mechanism = strings(kShow_astar,1);

    for i = 1:kShow_astar
        a = NodeA_astar(i);
        b = NodeB_astar(i);

        Dose_Free(i) = doseMat_free(a,b);
        Dose_Obs(i)  = doseMat_obs(a,b);

        Len_Free(i) = lenMat_free(a,b);
        Len_Obs(i)  = lenMat_obs(a,b);

        DoseChangePct(i) = 100 * (Dose_Obs(i) - Dose_Free(i)) / max(Dose_Free(i), eps);
        LenChangePct(i)  = 100 * (Len_Obs(i)  - Len_Free(i))  / max(Len_Free(i),  eps);

        if Dose_Obs(i) < Dose_Free(i) && Len_Obs(i) > Len_Free(i)
            Mechanism(i) = "增距降剂量";
        elseif Dose_Obs(i) > Dose_Free(i) && Len_Obs(i) > Len_Free(i)
            Mechanism(i) = "增距增剂量(主要避障)";
        elseif abs(Dose_Obs(i)-Dose_Free(i)) <= 1e-10 && Len_Obs(i) > Len_Free(i)
            Mechanism(i) = "增距剂量近似不变";
        else
            Mechanism(i) = "其他";
        end
    end

    Tastar = table(EdgeIndex_astar, NodeA_astar, NodeB_astar, ...
        Dose_Free, Dose_Obs, Len_Free, Len_Obs, ...
        DoseChangePct, LenChangePct, Mechanism);

    writetable(Tastar, fullfile(outDir, 'critical_astar_edge_compare.csv'));
end

%% ===================== 7) 关键变化边分析 =====================
edges_free = sort(tourToEdges(tour_free), 2);
edges_obs  = sort(tourToEdges(tour_obs),  2);

[onlyFree, onlyObs, ~] = compareEdgeSets(edges_free, edges_obs); %#ok<ASGLU>

if isempty(onlyObs)
    warning('没有发现障碍后新增边。');
    Tchange = table();
else
    nChange = size(onlyObs,1);

    NodeA_change = onlyObs(:,1);
    NodeB_change = onlyObs(:,2);

    Dose_Free_change = zeros(nChange,1);
    Dose_Obs_change  = zeros(nChange,1);
    Len_Free_change  = zeros(nChange,1);
    Len_Obs_change   = zeros(nChange,1);
    Type_Obs_change  = zeros(nChange,1);

    DoseDelta_change = zeros(nChange,1);
    LenDelta_change  = zeros(nChange,1);

    for i = 1:nChange
        a = NodeA_change(i);
        b = NodeB_change(i);

        Dose_Free_change(i) = doseMat_free(a,b);
        Dose_Obs_change(i)  = doseMat_obs(a,b);

        Len_Free_change(i) = lenMat_free(a,b);
        Len_Obs_change(i)  = lenMat_obs(a,b);

        Type_Obs_change(i) = edgeType_obs(a,b);

        DoseDelta_change(i) = Dose_Obs_change(i) - Dose_Free_change(i);
        LenDelta_change(i)  = Len_Obs_change(i)  - Len_Free_change(i);
    end

    Tchange = table(NodeA_change, NodeB_change, Type_Obs_change, ...
        Dose_Free_change, Dose_Obs_change, Len_Free_change, Len_Obs_change, ...
        DoseDelta_change, LenDelta_change);

    [~, idxSortChange] = sort(Dose_Obs_change, 'descend');
    Tchange = Tchange(idxSortChange,:);

    writetable(Tchange, fullfile(outDir, 'changed_edge_compare.csv'));
end

%% ===================== 8) 图1：总体高贡献边排序图 =====================
fig1 = figure('Color','w','Name','总体高贡献边排序图（论文风格）','Position',[110,110,1060,620]);
hold on; box on; grid on;

for i = 1:kShow_all
    if typeTop_all(i) == 0
        barh(i, doseTop_all(i), 0.62, 'FaceColor', col_los, 'EdgeColor', 'none');
    else
        barh(i, doseTop_all(i), 0.62, 'FaceColor', col_astar, 'EdgeColor', 'none');
    end
end

set(gca, 'YTick', 1:kShow_all, ...
    'YTickLabel', labelsTop_all, ...
    'YDir', 'reverse', ...
    'FontName', fontName, 'FontSize', fsTick, 'LineWidth', 1.1);

xlabel('剂量贡献 (a.u.)','FontName',fontName,'FontSize',fsAxis);
ylabel('边对（按剂量贡献排序）','FontName',fontName,'FontSize',fsAxis);
title(sprintf('有障碍混合模型中高剂量贡献边排序（Top-%d）', kShow_all), ...
    'FontName',fontName,'FontSize',fsTitle,'FontWeight','bold');

xMaxText = max(doseTop_all) * 1.02;
for i = 1:kShow_all
    if typeTop_all(i) == 0
        txt = sprintf('LOS | %.1f m', lenTop_all(i));
    else
        txt = sprintf('A* | %.1f m', lenTop_all(i));
    end
    text(xMaxText, i, txt, 'FontName',fontName, 'FontSize',10, ...
        'VerticalAlignment','middle');
end

legend({'LOS边','A*边'}, 'Location','southeast', ...
    'FontName',fontName,'FontSize',fsTick);

hold off;

saveas(fig1, fullfile(outDir, 'dose_breakdown_sorted_barh_paper.png'));

%% ===================== 9) 图2：关键A*边直线 vs 绕障对比散点 =====================
if ~isempty(idxAstar)
    fig2 = figure('Color','w','Name','关键A*边对比图（论文风格）','Position',[130,130,960,500]);
    hold on; box on; grid on;

    scatter(Len_Free, Dose_Free, 65, 'o', ...
        'MarkerFaceColor', col_free, 'MarkerEdgeColor', 'k', 'LineWidth', 0.8);
    scatter(Len_Obs, Dose_Obs, 75, 's', ...
        'MarkerFaceColor', col_astar, 'MarkerEdgeColor', 'k', 'LineWidth', 0.8);

    for i = 1:numel(Len_Free)
        plot([Len_Free(i), Len_Obs(i)], [Dose_Free(i), Dose_Obs(i)], ...
            '-', 'Color', [0.55 0.55 0.55], 'LineWidth', 1.2);
        text(Len_Obs(i), Dose_Obs(i), sprintf(' (%d,%d)', NodeA_astar(i), NodeB_astar(i)), ...
            'FontName', fontName, 'FontSize', 10);
    end

    xlabel('边长度 (m)','FontName',fontName,'FontSize',fsAxis);
    ylabel('边暴露代价 (a.u.)','FontName',fontName,'FontSize',fsAxis);
    title('关键 A* 边在无障碍直线与有障碍绕障条件下的长度-剂量变化', ...
        'FontName',fontName,'FontSize',fsTitle,'FontWeight','bold');
    legend({'无障碍直线边','有障碍A*边'}, 'Location','northwest', ...
        'FontName',fontName,'FontSize',fsTick);

    hold off;

    saveas(fig2, fullfile(outDir, 'critical_astar_compare_scatter_paper.png'));
end

%% ===================== 10) 图3：关键变化边对比条形图 =====================
if ~isempty(Tchange)
    kShow_change = min(topK_change, height(Tchange));
    TchangeTop = Tchange(1:kShow_change,:);

    labelsChange = strings(kShow_change,1);
    for i = 1:kShow_change
        labelsChange(i) = sprintf('(%d,%d)', TchangeTop.NodeA_change(i), TchangeTop.NodeB_change(i));
    end

    fig3 = figure('Color','w','Name','关键变化边对比图（论文风格）','Position',[140,140,1160,460]);

    subplot(1,2,1);
    hold on; box on; grid on;
    X = 1:kShow_change;
    bw = 0.35;
    bar(X-bw/2, TchangeTop.Dose_Free_change, bw, 'FaceColor', col_free, 'EdgeColor', 'none');
    bar(X+bw/2, TchangeTop.Dose_Obs_change,  bw, 'FaceColor', col_astar, 'EdgeColor', 'none');
    set(gca, 'XTick', X, 'XTickLabel', labelsChange, ...
        'FontName',fontName,'FontSize',fsTick,'LineWidth',1.1);
    xtickangle(35);
    ylabel('边暴露代价 (a.u.)','FontName',fontName,'FontSize',fsAxis);
    title('障碍后新增关键边：剂量比较', ...
        'FontName',fontName,'FontSize',fsTitle-1,'FontWeight','bold');
    legend({'无障碍直线','有障碍边'}, 'Location','best', ...
        'FontName',fontName,'FontSize',fsTick);
    hold off;

    subplot(1,2,2);
    hold on; box on; grid on;
    bar(X-bw/2, TchangeTop.Len_Free_change, bw, 'FaceColor', col_free, 'EdgeColor', 'none');
    bar(X+bw/2, TchangeTop.Len_Obs_change,  bw, 'FaceColor', col_astar, 'EdgeColor', 'none');
    set(gca, 'XTick', X, 'XTickLabel', labelsChange, ...
        'FontName',fontName,'FontSize',fsTick,'LineWidth',1.1);
    xtickangle(35);
    ylabel('边长度 (m)','FontName',fontName,'FontSize',fsAxis);
    title('障碍后新增关键边：长度比较', ...
        'FontName',fontName,'FontSize',fsTitle-1,'FontWeight','bold');
    legend({'无障碍直线','有障碍边'}, 'Location','best', ...
        'FontName',fontName,'FontSize',fsTick);
    hold off;

    sgtitle('障碍引入后关键替换边的长度与剂量变化', ...
        'FontName',fontName,'FontSize',fsTitle,'FontWeight','bold');

    saveas(fig3, fullfile(outDir, 'changed_edge_compare_bar_paper.png'));
end

%% ===================== 11) 关键A*边局部机制图：单独成图 =====================
if ~isempty(idxAstar)
    kLocal = min(height(Tastar), topK_astar);

    visRes = 2.0;
    xGrid = xMin:visRes:xMax;
    yGrid = yMin:visRes:yMax;
    [XX, YY] = meshgrid(xGrid, yGrid);
    radVis = zeros(size(XX));
    for s = 1:size(sources,1)
        dx = XX - sources(s,1);
        dy = YY - sources(s,2);
        radVis = radVis + sourceStrength(s) ./ (dx.^2 + dy.^2 + 10);
    end

    for i = 1:kLocal
        a = Tastar.NodeA_astar(i);
        b = Tastar.NodeB_astar(i);

        figLocal = figure('Color','w', ...
            'Name', sprintf('关键A*边局部机制图_(%d,%d)', a, b), ...
            'Position', [120,120,880,720]);
        hold on; box on; grid on; axis equal;

        contourf(XX, YY, radVis, 26, 'LineStyle','none');
        colormap(jet);
        cb = colorbar;
        cb.Label.String = '辐射水平（相对值）';
        cb.Label.FontName = fontName;
        cb.FontName = fontName;

        % 障碍物并标号
        drawObstacles_CN_withLabel(obstacles, col_obsFill, fontName);

        % 所有节点淡化显示
        scatter(xCoords, yCoords, 28, 'k', 'filled', 'MarkerEdgeColor','w', 'LineWidth',0.5);

        % 全部节点编号
        for n = 1:numel(xCoords)
            text(xCoords(n)+1.5, yCoords(n)+1.0, sprintf('%d', n), ...
                'Color','w', 'FontWeight','bold', ...
                'BackgroundColor',[0 0 0 0.30], ...
                'FontName',fontName, 'FontSize',8, ...
                'HorizontalAlignment','center');
        end

        % 点源
        scatter(sources(:,1), sources(:,2), 90, 'r', 'pentagram', 'filled');
        for s = 1:size(sources,1)
            text(sources(s,1), sources(s,2)-4, sprintf('源%d',s), ...
                'Color','r','FontWeight','bold', ...
                'HorizontalAlignment','center', ...
                'FontName',fontName,'FontSize',9);
        end

        % 无障碍直线
        plot([xCoords(a), xCoords(b)], [yCoords(a), yCoords(b)], ...
            '-', 'Color', col_free, 'LineWidth', 2.2);

        % A*真实绕障路径
        [pathRC, ~, ~] = astar_grid_dose_compare(nodeRC(a,:), nodeRC(b,:), ...
            occ, radGrid, gridRes, epsLen, minRad);
        pathRC2 = simplifyCollinear_compare(pathRC);
        pathRC3 = shortcutLineOfSight_compare(pathRC2, occ);
        [px, py] = gridPathToWorld_compare(pathRC3, xVec, yVec);

        % ===== 修正：首尾对齐到真实节点坐标 =====
        px(1)   = xCoords(a);
        py(1)   = yCoords(a);
        px(end) = xCoords(b);
        py(end) = yCoords(b);

        plot(px, py, '--', 'Color', col_astar, 'LineWidth', 2.5);

        % 高亮关键端点
        scatter([xCoords(a), xCoords(b)], [yCoords(a), yCoords(b)], ...
            90, 'y', 'filled', 'MarkerEdgeColor','k', 'LineWidth',1.0);

        text(xCoords(a), yCoords(a)+4, sprintf('节点%d', a), ...
            'Color','k', 'FontWeight','bold', ...
            'BackgroundColor','y', ...
            'HorizontalAlignment','center', ...
            'FontName',fontName,'FontSize',10);

        text(xCoords(b), yCoords(b)+4, sprintf('节点%d', b), ...
            'Color','k', 'FontWeight','bold', ...
            'BackgroundColor','y', ...
            'HorizontalAlignment','center', ...
            'FontName',fontName,'FontSize',10);

        title(sprintf('关键 A* 边局部路径机制图：边 (%d,%d)\n无障碍直线 vs 有障碍绕障 | 长度变化 %.1f%%，剂量变化 %.1f%%', ...
            a, b, Tastar.LenChangePct(i), Tastar.DoseChangePct(i)), ...
            'FontName',fontName,'FontSize',13,'FontWeight','bold');

        xlabel('X坐标','FontName',fontName,'FontSize',fsAxis);
        ylabel('Y坐标','FontName',fontName,'FontSize',fsAxis);
        set(gca,'FontName',fontName,'FontSize',fsTick,'LineWidth',1.1);

        % 局部范围
        xmid = (xCoords(a)+xCoords(b))/2;
        ymid = (yCoords(a)+yCoords(b))/2;
        xlim([max(xMin, xmid-localPad), min(xMax, xmid+localPad)]);
        ylim([max(yMin, ymid-localPad), min(yMax, ymid+localPad)]);

        % 图例
        h1 = plot(nan,nan,'-','Color',col_free,'LineWidth',2.2);
        h2 = plot(nan,nan,'--','Color',col_astar,'LineWidth',2.5);
        h3 = plot(nan,nan,'s','MarkerFaceColor',col_obsFill,'MarkerEdgeColor','k','LineStyle','none');
        h4 = plot(nan,nan,'o','MarkerFaceColor','y','MarkerEdgeColor','k','LineStyle','none');
        legend([h1,h2,h3,h4], {'无障碍直线边','有障碍A*边','障碍物','关键端点'}, ...
            'Location','best', 'FontName',fontName, 'FontSize',10);

        hold off;

        saveas(figLocal, fullfile(outDir, sprintf('local_key_edge_%02d_(%d_%d).png', i, a, b)));
    end
end

%% ===================== 12) 控制台总结 =====================
fprintf('\n===== 剂量分解分析结果 =====\n');
fprintf('有障碍最优路径中前 %d 条高剂量贡献边已输出。\n', kShow_all);

if ~isempty(idxAstar)
    fprintf('关键 A* 边分析数量: %d\n', height(Tastar));
    mechSummary = groupsummary(Tastar, 'Mechanism');
    disp(mechSummary);
else
    fprintf('有障碍最优路径中无 A* 边。\n');
end

if ~isempty(Tchange)
    fprintf('障碍后新增关键边数量: %d\n', height(Tchange));
end

%% ===================== 13) 保存 MAT =====================
save(fullfile(outDir, 'analyze_dose_breakdown_case3_paper_v2.mat'), ...
    'Tall','TtopAll','Tastar','Tchange', ...
    'doseList_obs','lenList_obs','typeList_obs','edgePairs_obs', ...
    'idxSort_all','idxTop_all','doseTop_all','lenTop_all','typeTop_all','pairTop_all');

fprintf('\n✅ 剂量分解分析脚本执行完成。\n');
fprintf('输出目录：%s\n', outDir);

end


%% ============================================================
%                          Functions
%% ============================================================

function [doseList, lenList, typeList, edgePairs] = getPathEdgeDetail(tour, doseMat, lenMat, edgeType)
m = numel(tour)-1;
doseList = zeros(m,1);
lenList  = zeros(m,1);
typeList = zeros(m,1);
edgePairs = zeros(m,2);

for i = 1:m
    a = tour(i);
    b = tour(i+1);
    doseList(i) = doseMat(a,b);
    lenList(i)  = lenMat(a,b);
    typeList(i) = edgeType(a,b);
    edgePairs(i,:) = [a,b];
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

function drawObstacles_CN_withLabel(obstacles, faceColor, fontName)
for k = 1:size(obstacles,1)
    r = obstacles(k,:);
    rectangle('Position', r, 'FaceColor', [faceColor, 0.65], ...
        'EdgeColor','k', 'LineWidth',1.2);
    cx = r(1) + r(3)/2;
    cy = r(2) + r(4)/2;
    text(cx, cy, sprintf('障碍%d', k), ...
        'Color','w','FontWeight','bold', ...
        'HorizontalAlignment','center', ...
        'VerticalAlignment','middle', ...
        'FontName',fontName,'FontSize',10);
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
    error('A* failed in dose breakdown script: (%d,%d)->(%d,%d) unreachable.', sr,sc,gr,gc);
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
        error('Path reconstruction failed in dose breakdown script.');
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