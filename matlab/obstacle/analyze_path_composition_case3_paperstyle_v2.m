function analyze_path_composition_case3_paperstyle_v2()
% ============================================================
% analyze_path_composition_case3_paperstyle_v2.m
% ------------------------------------------------------------
% 第三部分分析脚本 2（升级版）：
% 有障碍混合模型下的路径组成分析（论文风格增强）
%
% 输入：
%   1) result_case3_LOS_refined.mat
%   2) result_case3_mixed_obstacle.mat
%
% 输出：
%   1) path_composition_summary.csv
%   2) path_edge_detail_mixed.csv
%   3) composition_absolute_bar_paper.png
%   4) composition_ratio_stacked_paper.png
%   5) top_contributing_edges_barh_paper.png
%   6) analyze_path_composition_case3_paper_v2.mat
%
% 核心分析：
%   - LOS / A* 边数量
%   - LOS / A* 剂量总和
%   - LOS / A* 距离总和
%   - A*边剂量占比 / 距离占比
%   - 高贡献边排序分析（论文增强版）
% ============================================================

clc; close all;

%% ===================== 0) 文件路径 =====================
fileLOS   = fullfile(pwd, 'result_case3_LOS_refined.mat');
fileMixed = fullfile(pwd, 'result_case3_mixed_obstacle.mat');

assert(exist(fileLOS, 'file')==2, '找不到文件: %s', fileLOS);
assert(exist(fileMixed, 'file')==2, '找不到文件: %s', fileMixed);

outDir = fullfile(pwd, 'analyze_path_composition_case3_paper_v2');
if ~exist(outDir, 'dir')
    mkdir(outDir);
end

fprintf('\n================= 路径组成分析（升级版，论文风格） =================\n');
fprintf('LOS file   : %s\n', fileLOS);
fprintf('Mixed file : %s\n', fileMixed);
fprintf('Output dir : %s\n', outDir);

%% ===================== 1) 读取数据 =====================
A = load(fileLOS);
B = load(fileMixed);

reqLOS = {'gBestPath','doseMat','lenMat','edgeType','totalDose','totalLen'};
reqMix = {'gBestPath','doseMat','lenMat','edgeType','totalDose','totalLen'};

for i = 1:numel(reqLOS)
    assert(isfield(A, reqLOS{i}), '无障碍结果缺字段: %s', reqLOS{i});
end
for i = 1:numel(reqMix)
    assert(isfield(B, reqMix{i}), '有障碍结果缺字段: %s', reqMix{i});
end

%% ===================== 2) 论文风格设置 =====================
fontName = 'SimSun';
fsAxis   = 12;
fsTitle  = 14;
fsTick   = 11;

col_los   = [0.10, 0.35, 0.80];
col_astar = [0.80, 0.15, 0.15];
col_gray  = [0.45, 0.45, 0.45];
col_free  = [0.15, 0.55, 0.55];

topK = 15;   % 高贡献边显示前K条

%% ===================== 3) 提取最优路径信息 =====================
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

%% ===================== 4) 逐边细节 =====================
[doseList_free, lenList_free, typeList_free, edgePairs_free] = ...
    getPathCompositionDetail(tour_free, doseMat_free, lenMat_free, edgeType_free);

[doseList_obs, lenList_obs, typeList_obs, edgePairs_obs] = ...
    getPathCompositionDetail(tour_obs, doseMat_obs, lenMat_obs, edgeType_obs);

%% ===================== 5) 无障碍模型组成 =====================
nLOS_free   = sum(typeList_free == 0);
nAStar_free = sum(typeList_free == 1);

doseLOS_free   = sum(doseList_free(typeList_free == 0));
doseAStar_free = sum(doseList_free(typeList_free == 1));

lenLOS_free   = sum(lenList_free(typeList_free == 0));
lenAStar_free = sum(lenList_free(typeList_free == 1));

ratioDoseAStar_free = 100 * doseAStar_free / max(totalDose_free, eps);
ratioLenAStar_free  = 100 * lenAStar_free  / max(totalLen_free,  eps);

%% ===================== 6) 有障碍模型组成 =====================
nLOS_obs   = sum(typeList_obs == 0);
nAStar_obs = sum(typeList_obs == 1);

doseLOS_obs   = sum(doseList_obs(typeList_obs == 0));
doseAStar_obs = sum(doseList_obs(typeList_obs == 1));

lenLOS_obs   = sum(lenList_obs(typeList_obs == 0));
lenAStar_obs = sum(lenList_obs(typeList_obs == 1));

ratioDoseAStar_obs = 100 * doseAStar_obs / max(totalDose_obs, eps);
ratioLenAStar_obs  = 100 * lenAStar_obs  / max(totalLen_obs,  eps);

ratioDoseLOS_obs = 100 - ratioDoseAStar_obs;
ratioLenLOS_obs  = 100 - ratioLenAStar_obs;

%% ===================== 7) 导出汇总表 =====================
Model = ["无障碍精细模型"; "有障碍混合模型"];

LOSEdgeCount   = [nLOS_free; nLOS_obs];
AStarEdgeCount = [nAStar_free; nAStar_obs];

LOSEdgeDose   = [doseLOS_free; doseLOS_obs];
AStarEdgeDose = [doseAStar_free; doseAStar_obs];

LOSEdgeLength   = [lenLOS_free; lenLOS_obs];
AStarEdgeLength = [lenAStar_free; lenAStar_obs];

AStarDoseRatioPct   = [ratioDoseAStar_free; ratioDoseAStar_obs];
AStarLengthRatioPct = [ratioLenAStar_free; ratioLenAStar_obs];

Tsummary = table(Model, ...
    LOSEdgeCount, AStarEdgeCount, ...
    LOSEdgeDose, AStarEdgeDose, ...
    LOSEdgeLength, AStarEdgeLength, ...
    AStarDoseRatioPct, AStarLengthRatioPct);

writetable(Tsummary, fullfile(outDir, 'path_composition_summary.csv'));

fprintf('\n===== 路径组成汇总 =====\n');
disp(Tsummary);

%% ===================== 8) 导出有障碍逐边明细表 =====================
EdgeIndex = (1:numel(doseList_obs)).';
NodeA = edgePairs_obs(:,1);
NodeB = edgePairs_obs(:,2);

EdgeType = strings(numel(doseList_obs),1);
EdgeType(typeList_obs==0) = "LOS";
EdgeType(typeList_obs==1) = "A*";

DoseRatioPct = 100 * doseList_obs / max(totalDose_obs, eps);
LengthRatioPct = 100 * lenList_obs / max(totalLen_obs, eps);

Tdetail = table(EdgeIndex, NodeA, NodeB, EdgeType, doseList_obs, lenList_obs, ...
    DoseRatioPct, LengthRatioPct, ...
    'VariableNames', {'EdgeIndex','NodeA','NodeB','EdgeType','Dose','Length','DoseRatioPct','LengthRatioPct'});

writetable(Tdetail, fullfile(outDir, 'path_edge_detail_mixed.csv'));

%% ===================== 9) 图1：绝对量组成柱状图 =====================
fig1 = figure('Color','w','Name','路径组成绝对量柱状图（论文风格）','Position',[120,120,1150,440]);

subplot(1,3,1);
hold on; box on; grid on;
vals_count = [nLOS_obs, nAStar_obs];
b = bar(vals_count, 0.58);
b.FaceColor = 'flat';
b.CData = [col_los; col_astar];
set(gca,'XTick',1:2,'XTickLabel',{'LOS边','A*边'}, ...
    'FontName',fontName,'FontSize',fsTick,'LineWidth',1.1);
ylabel('边数量','FontName',fontName,'FontSize',fsAxis);
title('边数量组成','FontName',fontName,'FontSize',fsTitle-1,'FontWeight','bold');
hold off;

subplot(1,3,2);
hold on; box on; grid on;
vals_dose = [doseLOS_obs, doseAStar_obs];
b = bar(vals_dose, 0.58);
b.FaceColor = 'flat';
b.CData = [col_los; col_astar];
set(gca,'XTick',1:2,'XTickLabel',{'LOS边','A*边'}, ...
    'FontName',fontName,'FontSize',fsTick,'LineWidth',1.1);
ylabel('剂量贡献 (\muSv)','FontName',fontName,'FontSize',fsAxis);
title('剂量绝对贡献','FontName',fontName,'FontSize',fsTitle-1,'FontWeight','bold');
hold off;

subplot(1,3,3);
hold on; box on; grid on;
vals_len = [lenLOS_obs, lenAStar_obs];
b = bar(vals_len, 0.58);
b.FaceColor = 'flat';
b.CData = [col_los; col_astar];
set(gca,'XTick',1:2,'XTickLabel',{'LOS边','A*边'}, ...
    'FontName',fontName,'FontSize',fsTick,'LineWidth',1.1);
ylabel('距离贡献 (m)','FontName',fontName,'FontSize',fsAxis);
title('距离绝对贡献','FontName',fontName,'FontSize',fsTitle-1,'FontWeight','bold');
hold off;

sgtitle('有障碍混合模型最优路径的绝对量组成分析', ...
    'FontName',fontName,'FontSize',fsTitle,'FontWeight','bold');

saveas(fig1, fullfile(outDir, 'composition_absolute_bar_paper.png'));

%% ===================== 10) 图2：100%%堆叠柱状图（替代饼图） =====================
fig2 = figure('Color','w','Name','组成占比堆叠图（论文风格）','Position',[140,140,900,460]);
hold on; box on; grid on;

Y = [
    ratioDoseLOS_obs, ratioDoseAStar_obs;
    ratioLenLOS_obs,  ratioLenAStar_obs
];

b = bar(Y, 'stacked', 'BarWidth', 0.58);
b(1).FaceColor = col_los;
b(2).FaceColor = col_astar;

set(gca, 'XTick', 1:2, ...
    'XTickLabel', {'剂量占比', '距离占比'}, ...
    'FontName', fontName, 'FontSize', fsTick, 'LineWidth', 1.1);
ylabel('占比 (%)','FontName',fontName,'FontSize',fsAxis);
title('有障碍混合模型中 LOS / A* 边的占比分析', ...
    'FontName',fontName,'FontSize',fsTitle,'FontWeight','bold');
legend({'LOS边','A*边'}, 'Location','northeast', ...
    'FontName',fontName,'FontSize',fsTick);

% 标注比例
for i = 1:size(Y,1)
    y1 = Y(i,1)/2;
    y2 = Y(i,1) + Y(i,2)/2;
    text(i, y1, sprintf('%.1f%%', Y(i,1)), ...
        'HorizontalAlignment','center','Color','w', ...
        'FontName',fontName,'FontSize',10,'FontWeight','bold');
    text(i, y2, sprintf('%.1f%%', Y(i,2)), ...
        'HorizontalAlignment','center','Color','w', ...
        'FontName',fontName,'FontSize',10,'FontWeight','bold');
end

hold off;

saveas(fig2, fullfile(outDir, 'composition_ratio_stacked_paper.png'));

%% ===================== 11) 图3：高贡献边排序横向条形图 =====================
% 按剂量贡献从高到低排序
[sortedDose, idxSort] = sort(doseList_obs, 'descend');
kShow = min(topK, numel(sortedDose));

idxTop = idxSort(1:kShow);
doseTop = doseList_obs(idxTop);
lenTop  = lenList_obs(idxTop);
typeTop = typeList_obs(idxTop);
pairTop = edgePairs_obs(idxTop,:);

labelsTop = strings(kShow,1);
for i = 1:kShow
    labelsTop(i) = sprintf('(%d,%d)', pairTop(i,1), pairTop(i,2));
end

fig3 = figure('Color','w','Name','高贡献边排序图（论文风格）','Position',[120,120,1050,620]);
hold on; box on; grid on;

for i = 1:kShow
    if typeTop(i) == 0
        barh(i, doseTop(i), 0.62, 'FaceColor', col_los, 'EdgeColor', 'none');
    else
        barh(i, doseTop(i), 0.62, 'FaceColor', col_astar, 'EdgeColor', 'none');
    end
end

set(gca, 'YTick', 1:kShow, ...
    'YTickLabel', labelsTop, ...
    'YDir', 'reverse', ...
    'FontName', fontName, 'FontSize', fsTick, 'LineWidth', 1.1);

xlabel('剂量贡献 (\muSv)','FontName',fontName,'FontSize',fsAxis);
ylabel('边对（按剂量贡献排序）','FontName',fontName,'FontSize',fsAxis);
title(sprintf('有障碍混合模型中高贡献边排序分析（Top-%d）', kShow), ...
    'FontName',fontName,'FontSize',fsTitle,'FontWeight','bold');

% 右侧标注类型和长度
xMaxText = max(doseTop) * 1.02;
for i = 1:kShow
    if typeTop(i) == 0
        tStr = sprintf('LOS | %.1f m', lenTop(i));
    else
        tStr = sprintf('A* | %.1f m', lenTop(i));
    end
    text(xMaxText, i, tStr, 'FontName',fontName, 'FontSize',12, ...
        'VerticalAlignment','middle');
end

legend({'LOS边','A*边'}, 'Location','southeast', ...
    'FontName',fontName,'FontSize',fsTick);

hold off;

saveas(fig3, fullfile(outDir, 'top_contributing_edges_barh_paper.png'));

%% ===================== 12) 输出文字总结 =====================
fprintf('\n===== 有障碍路径组成核心结果 =====\n');
fprintf('LOS边数量: %d\n', nLOS_obs);
fprintf('A*边数量: %d\n', nAStar_obs);
fprintf('LOS边剂量总和: %.4f\n', doseLOS_obs);
fprintf('A*边剂量总和: %.4f\n', doseAStar_obs);
fprintf('LOS边距离总和: %.4f\n', lenLOS_obs);
fprintf('A*边距离总和: %.4f\n', lenAStar_obs);
fprintf('A*边剂量占总剂量比例: %.2f%%\n', ratioDoseAStar_obs);
fprintf('A*边距离占总距离比例: %.2f%%\n', ratioLenAStar_obs);

%% ===================== 13) 保存 MAT =====================
save(fullfile(outDir, 'analyze_path_composition_case3_paper_v2.mat'), ...
    'tour_free','tour_obs', ...
    'doseList_free','lenList_free','typeList_free','edgePairs_free', ...
    'doseList_obs','lenList_obs','typeList_obs','edgePairs_obs', ...
    'nLOS_free','nAStar_free','doseLOS_free','doseAStar_free','lenLOS_free','lenAStar_free', ...
    'nLOS_obs','nAStar_obs','doseLOS_obs','doseAStar_obs','lenLOS_obs','lenAStar_obs', ...
    'ratioDoseAStar_free','ratioLenAStar_free', ...
    'ratioDoseAStar_obs','ratioLenAStar_obs', ...
    'ratioDoseLOS_obs','ratioLenLOS_obs', ...
    'Tsummary','Tdetail', ...
    'idxSort','idxTop','doseTop','lenTop','typeTop','pairTop');

fprintf('\n✅ 升级版路径组成分析脚本执行完成。\n');
fprintf('输出目录：%s\n', outDir);

end


%% ============================================================
%                          Functions
%% ============================================================

function [doseList, lenList, typeList, edgePairs] = getPathCompositionDetail(tour, doseMat, lenMat, edgeType)
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
    edgePairs(i,:) = [a, b];
end
end