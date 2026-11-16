function benchmark_merge_all_algorithms_paperstyle_final_CN_Bcolors()
% benchmark_merge_all_algorithms_paperstyle_final_CN_Bcolors.m
% ============================================================
%  总控汇总（中文增强版 + 方案B散点按算法着色 + 不自动保存png）
%  读取各算法 results_*_20instances.mat，生成：
%
%  (A) CSV（输出到 ./bench_ALL_merged_paper_final_CN/ ）：
%      1) ranked_results_ALL_algorithms.csv  (长表：Algorithm-Instance一行，按BestDose排序)
%      2) summary_stats_ALL_algorithms.csv   (每算法：mean/std/median/p90/p95/max 等)
%      3) pairwise_mannwhitney_BestDose.csv  (两两Mann-Whitney + BH校正 + A12效应量)
%      4) text_summary_ALL_algorithms.txt    (新增：文字总结)
%
%  (B) 图（窗口形式打开，你手动调整后保存）：
%      1) 多算法并排箱线图+散点（仅普通标记；散点按算法着色=方案B）：
%         - BestDose / BestLen / Time
%      2) 多算法均值±1σ 收敛对比图（含两个放大窗）
%      3) ECDF（BestDose）多算法对比
%      4) FirstHitIter / LastImproveIter 并排箱线图+散点（散点按算法着色）
%      5) 速度-质量权衡散点：x=FirstHitIter, y=BestDose
%      6) 噪声敏感性：BestDose vs Noise / Runtime vs Noise（若有Noise）
%      7) 最差情形风险条形图：P90/P95/Max
%      8) （可选）综合评分条形图（Dose/Time/Len 归一化加权）
% ============================================================

clc; close all;

%% ===================== 0) 配置：算法与文件路径 =====================
root = pwd;

algoList = { ...
    struct('name',"CPO",     'mat', fullfile(root,'bench_CPO_20_paper','results_CPO_20instances.mat')), ...
    struct('name',"DBO",     'mat', fullfile(root,'bench_DBO_20_paper','results_DBO_20instances.mat')), ...
    struct('name',"CPO-DBO", 'mat', fullfile(root,'bench_CPO_DBO_20_paper','results_CPO_DBO_20instances.mat')), ...
    struct('name',"IPSO",    'mat', fullfile(root,'bench_IPSO_20_paper','results_IPSO_20instances.mat')), ...
    struct('name',"GA",      'mat', fullfile(root,'bench_GA_20_paper','results_GA_20instances.mat')) ...
};

outDir = fullfile(root, 'bench_ALL_merged_paper_final_CN');
if ~exist(outDir,'dir'), mkdir(outDir); end

theme = "orange";         % "orange" / "purpleblue"
doCompositeScore = true;  % 是否生成综合评分图

% 综合评分权重（可改）
wDose = 0.5;
wTime = 0.1;
wLen  = 0.4;

fprintf('\n================= 总控汇总（中文增强版 + 方案B散点着色）=================\n');
fprintf('CSV输出目录：%s\n', outDir);

%% ===================== 1) 读取并校验每个算法结果 =====================
A = numel(algoList);
All = struct();

for a = 1:A
    f = algoList{a}.mat;
    assert(exist(f,'file')==2, '找不到结果文件：%s', f);
    S = load(f);

    must = {'bestDoseRuns','bestLenRuns','timeRuns','convAll','maxIter','numInstances','firstHitIter','lastImproveIter'};
    for k = 1:numel(must)
        assert(isfield(S, must{k}), '文件 %s 缺字段：%s', f, must{k});
    end

    algoName = algoList{a}.name;
    All.(matlab.lang.makeValidName(algoName)) = S;

    fprintf('Loaded %-7s | instances=%d | maxIter=%d | %s\n', algoName, S.numInstances, S.maxIter, getFileName(f));
end

numInstances = All.(matlab.lang.makeValidName(algoList{1}.name)).numInstances;
for a = 2:A
    s = All.(matlab.lang.makeValidName(algoList{a}.name));
    assert(s.numInstances == numInstances, '算法 %s 的 numInstances 不一致', algoList{a}.name);
end

algoLabels = cellfun(@(s)char(s.name), algoList, 'UniformOutput', false);

%% ===================== 2) 生成“长表”论文总表（并导出CSV） =====================
rows = A * numInstances;
Algorithm  = strings(rows,1);
Instance   = zeros(rows,1);
BestDose   = zeros(rows,1);
BestLen    = zeros(rows,1);
FirstHit   = zeros(rows,1);
LastImp    = zeros(rows,1);
TimeSec    = zeros(rows,1);
NoiseLevel = nan(rows,1);
SceneFP    = nan(rows,1);

r = 0;
for a = 1:A
    algoName = algoList{a}.name;
    S = All.(matlab.lang.makeValidName(algoName));

    for k = 1:numInstances
        r = r + 1;
        Algorithm(r) = algoName;
        Instance(r)  = k;
        BestDose(r)  = S.bestDoseRuns(k);
        BestLen(r)   = S.bestLenRuns(k);
        FirstHit(r)  = S.firstHitIter(k);
        LastImp(r)   = S.lastImproveIter(k);
        TimeSec(r)   = S.timeRuns(k);

        if isfield(S,'instNoiseLevel'), NoiseLevel(r) = S.instNoiseLevel(k); end
        if isfield(S,'sceneFP'),        SceneFP(r)    = S.sceneFP(k);        end
    end
end

Tall = table(Algorithm, Instance, BestDose, BestLen, FirstHit, LastImp, TimeSec, NoiseLevel, SceneFP);
[~, idxSortAll] = sort(Tall.BestDose, 'ascend');
Tall = Tall(idxSortAll, :);
Tall.Rank = (1:height(Tall))';
Tall = movevars(Tall, 'Rank', 'Before', 1);

outCsv = fullfile(outDir, 'ranked_results_ALL_algorithms.csv');
writetable(Tall, outCsv);
fprintf('\n✅ 已导出论文总表：%s\n', outCsv);

%% ===================== 3) 论文风样式 =====================
switch theme
    case "orange"
        shadeColor = [0.95, 0.45, 0.10];
        meanColor  = [0.05, 0.05, 0.05];
    case "purpleblue"
        shadeColor = [0.35, 0.45, 0.95];
        meanColor  = [0.05, 0.05, 0.05];
end
shadeAlpha = 0.18;
lwMean = 2.2;
lwBox = 1.6;
lwMedian = 2.2;
jit = 0.12;
mkSize = 22;
mkLW = 1.2;

fontName = 'SimSun';
fsAxis = 12;
fsTitle = 14;

grayEdge = [0.25 0.25 0.25];
grayFace = [0.20 0.20 0.20];

%% ===================== 4) 图1：多算法并排箱线图+散点（方案B：按算法着色） =====================
makeMultiAlgoBoxScatter_plain_colorB(All, algoList, numInstances, 'bestDoseRuns', ...
    '最优总剂量（\muSv）', '多算法对比：最优总剂量（箱线图+散点）', ...
    algoLabels, fontName, fsAxis, fsTitle, lwBox, lwMedian, jit, mkSize, mkLW, grayEdge, grayFace);

makeMultiAlgoBoxScatter_plain_colorB(All, algoList, numInstances, 'bestLenRuns', ...
    '最优路径总距离（m）', '多算法对比：最优路径总距离（箱线图+散点）', ...
    algoLabels, fontName, fsAxis, fsTitle, lwBox, lwMedian, jit, mkSize, mkLW, grayEdge, grayFace);

makeMultiAlgoBoxScatter_plain_colorB(All, algoList, numInstances, 'timeRuns', ...
    '运行时间（s）', '多算法对比：运行时间（箱线图+散点）', ...
    algoLabels, fontName, fsAxis, fsTitle, lwBox, lwMedian, jit, mkSize, mkLW, grayEdge, grayFace);

%% ===================== 5) 图2：多算法均值±1σ 收敛对比（含两个放大窗） =====================
zoom1 = [1,   400,  320, 700];
zoom2 = [2500, 3000, 320, 480];

makeConvergenceMeanStdCompare_2inset(All, algoList, numInstances, ...
    fontName, fsAxis, fsTitle, lwMean, shadeAlpha, zoom1, zoom2);

%% ===================== 6) 图3：ECDF（BestDose）多算法对比 =====================
makeECDFCompare(All, algoList, numInstances, fontName, fsAxis, fsTitle);

%% ===================== 7) 图4：收敛速度指标箱线图+散点（方案B：按算法着色） =====================
makeMultiAlgoBoxScatter_plain_colorB(All, algoList, numInstances, 'firstHitIter', ...
    'FirstHitIter（迭代步）', '多算法对比：首次达到最终最优的迭代步（箱线图+散点）', ...
    algoLabels, fontName, fsAxis, fsTitle, lwBox, lwMedian, jit, mkSize, mkLW, grayEdge, grayFace);

makeMultiAlgoBoxScatter_plain_colorB(All, algoList, numInstances, 'lastImproveIter', ...
    'LastImproveIter（迭代步）', '多算法对比：最后一次改进的迭代步（箱线图+散点）', ...
    algoLabels, fontName, fsAxis, fsTitle, lwBox, lwMedian, jit, mkSize, mkLW, grayEdge, grayFace);

%% ===================== 8) 图5：速度-质量权衡散点（FirstHit vs BestDose） =====================
makeSpeedQualityScatter(All, algoList, numInstances, fontName, fsAxis, fsTitle);

%% ===================== 9) 图6：噪声敏感性（BestDose/Time vs Noise） =====================
hasAnyNoise = true;
for a = 1:A
    S = All.(matlab.lang.makeValidName(algoList{a}.name));
    if ~isfield(S,'instNoiseLevel') || all(isnan(S.instNoiseLevel(:)))
        hasAnyNoise = false;
    end
end

if hasAnyNoise
    makeNoiseScatterFit(All, algoList, numInstances, 'bestDoseRuns', ...
        '噪声水平（NoiseLevel）', 'BestDose（\muSv）', '噪声敏感性：BestDose 与 NoiseLevel 的关系', ...
        fontName, fsAxis, fsTitle);

    makeNoiseScatterFit(All, algoList, numInstances, 'timeRuns', ...
        '噪声水平（NoiseLevel）', '运行时间（s）', '噪声敏感性：运行时间 与 NoiseLevel 的关系', ...
        fontName, fsAxis, fsTitle);
else
    fprintf('\n[提示] 结果 MAT 中缺少 instNoiseLevel（或全为NaN），已跳过 Noise 敏感性图。\n');
end

%% ===================== 10) 图7：最差情形风险（P90/P95/Max） =====================
makeRiskBar(All, algoList, numInstances, fontName, fsAxis, fsTitle);

%% ===================== 11) （可选）综合评分条形图（归一化加权） =====================
if doCompositeScore
    makeCompositeScore(All, algoList, wDose, wTime, wLen, fontName, fsAxis, fsTitle);
end

%% ===================== 12) 导出每算法统计摘要 + 两两检验表 =====================
TsummaryStats = exportSummaryStats(All, algoList, outDir);
Tpairwise     = exportPairwiseStats_BestDose(All, algoList, outDir);

%% ===================== 13) 输出文字结果总结 =====================
printMergedTextSummary(All, algoList, numInstances, TsummaryStats, Tpairwise, ...
    outDir, doCompositeScore, wDose, wTime, wLen);

fprintf('\n===== 总控汇总（中文增强版 + 方案B散点着色）完成 =====\n');
fprintf('CSV输出目录：%s\n', outDir);
fprintf('图已全部以窗口形式生成（你可手动微调后保存）。\n');

end % main


%% ============================================================
%  图：多算法并排箱线图+散点（仅普通标记；方案B：按算法着色）
%% ============================================================
function makeMultiAlgoBoxScatter_plain_colorB(All, algoList, numInstances, fieldName, ylab, ttl, labels, ...
    fontName, fsAxis, fsTitle, lwBox, lwMedian, jit, mkSize, mkLW, grayEdge, grayFace)

A = numel(algoList);

data = cell(A,1);
for a = 1:A
    algoName = algoList{a}.name;
    S = All.(matlab.lang.makeValidName(algoName));
    vec = S.(fieldName)(:);
    assert(numel(vec)==numInstances, '字段 %s 的实例数不一致', fieldName);
    data{a} = vec;
end

Y = []; G = [];
for a = 1:A
    Y = [Y; data{a}]; %#ok<AGROW>
    G = [G; a*ones(numInstances,1)]; %#ok<AGROW>
end

figure('Color','w','Name',ttl,'Position',[120,120,1180,520]);
hold on; grid on; box on;

boxplot(Y, G, 'Labels', labels, 'Symbol','', 'Whisker',1.5);
set(findobj(gca,'Tag','Box'),     'LineWidth', lwBox);
set(findobj(gca,'Tag','Whisker'), 'LineWidth', lwBox);
set(findobj(gca,'Tag','Cap'),     'LineWidth', lwBox);
set(findobj(gca,'Tag','Median'),  'LineWidth', lwMedian);

hBox = findobj(gca,'Tag','Box');
for kk = 1:numel(hBox)
    xd = get(hBox(kk),'XData'); yd = get(hBox(kk),'YData');
    patch(xd, yd, grayFace, 'FaceAlpha',0.08, 'EdgeColor', grayEdge, 'LineWidth', lwBox);
end

cols = lines(A);

for a = 1:A
    vec = data{a};
    x0 = a;
    xj = x0 + (rand(numInstances,1)-0.5)*2*jit;

    scatter(xj, vec, mkSize, 'o', ...
        'MarkerFaceColor','none', ...
        'MarkerEdgeColor', cols(a,:), ...
        'LineWidth', mkLW);
end

ylabel(ylab,'FontName',fontName,'FontSize',fsAxis);
title(ttl, 'FontName',fontName,'FontSize',fsTitle,'FontWeight','bold');
set(gca,'FontName',fontName,'FontSize',fsAxis,'LineWidth',1.2);

hold off;
end


%% ============================================================
%  图：多算法均值±1σ 收敛对比（含两个放大窗）
%% ============================================================
function makeConvergenceMeanStdCompare_2inset(All, algoList, numInstances, ...
    fontName, fsAxis, fsTitle, lwMean, shadeAlpha, zoom1, zoom2)

A = numel(algoList);
cols = lines(A);

minIter = inf;
for a = 1:A
    S = All.(matlab.lang.makeValidName(algoList{a}.name));
    minIter = min(minIter, size(S.convAll,1));
end
iters = (1:minIter)';

muAll = zeros(minIter, A);
sdAll = zeros(minIter, A);

for a = 1:A
    algoName = algoList{a}.name;
    S = All.(matlab.lang.makeValidName(algoName));
    conv = S.convAll(1:minIter, :);
    muAll(:,a) = mean(conv, 2, 'omitnan');
    sdAll(:,a) = std(conv, 0, 2, 'omitnan');
end

fig = figure('Color','w','Name','多算法收敛对比（均值±1σ，含两个放大窗）','Position',[120,120,1220,560]);
axMain = axes(fig); %#ok<LAXES>
hold(axMain,'on'); grid(axMain,'on'); box(axMain,'on');

meanHandles = gobjects(A,1);
for a = 1:A
    xShade = [iters; flipud(iters)];
    yShade = [muAll(:,a)-sdAll(:,a); flipud(muAll(:,a)+sdAll(:,a))];
    fill(axMain, xShade, yShade, cols(a,:), 'EdgeColor','none', 'FaceAlpha', shadeAlpha);
    meanHandles(a) = plot(axMain, iters, muAll(:,a), '-', 'LineWidth', lwMean, 'Color', cols(a,:));
end

xlabel(axMain, '迭代次数','FontName',fontName,'FontSize',fsAxis);
ylabel(axMain, '最优剂量（\muSv）','FontName',fontName,'FontSize',fsAxis);
title(axMain, sprintf('多算法收敛对比（均值±1σ）：%d个扰动实例（对齐到 minIter=%d）', numInstances, minIter), ...
    'FontName',fontName,'FontSize',fsTitle,'FontWeight','bold');
set(axMain,'FontName',fontName,'FontSize',fsAxis,'LineWidth',1.2);

leg = cell(A,1);
for a = 1:A, leg{a} = char(algoList{a}.name); end
legend(axMain, meanHandles, leg, 'Location','northeastoutside');

zoom1 = resolveZoomY(zoom1, iters, muAll);
zoom2 = resolveZoomY(zoom2, iters, muAll);

rectangle(axMain, 'Position', [zoom1(1), zoom1(3), zoom1(2)-zoom1(1), zoom1(4)-zoom1(3)], ...
    'EdgeColor','r', 'LineWidth',1.8,'LineStyle','--');

rectangle(axMain, 'Position', [zoom2(1), zoom2(3), zoom2(2)-zoom2(1), zoom2(4)-zoom2(3)], ...
    'EdgeColor','r', 'LineWidth',1.8,'LineStyle','--');

insetPos1 = [0.2 0.42 0.18 0.22];
ax1 = axes('Position', insetPos1);
box(ax1,'on'); grid(ax1,'on'); hold(ax1,'on');
for a = 1:A
    xShade = [iters; flipud(iters)];
    yShade = [muAll(:,a)-sdAll(:,a); flipud(muAll(:,a)+sdAll(:,a))];
    fill(ax1, xShade, yShade, cols(a,:), 'EdgeColor','none', 'FaceAlpha', shadeAlpha);
    plot(ax1, iters, muAll(:,a), '-', 'LineWidth', 1.8, 'Color', cols(a,:));
end
xlim(ax1, [zoom1(1) zoom1(2)]);
ylim(ax1, [zoom1(3) zoom1(4)]);
set(ax1,'FontSize',10);
title(ax1,'放大区域1','FontSize',10,'FontWeight','bold');

insetPos2 = [0.58 0.42 0.18 0.28];
ax2 = axes('Position', insetPos2);
box(ax2,'on'); grid(ax2,'on'); hold(ax2,'on');
for a = 1:A
    xShade = [iters; flipud(iters)];
    yShade = [muAll(:,a)-sdAll(:,a); flipud(muAll(:,a)+sdAll(:,a))];
    fill(ax2, xShade, yShade, cols(a,:), 'EdgeColor','none', 'FaceAlpha', shadeAlpha);
    plot(ax2, iters, muAll(:,a), '-', 'LineWidth', 1.8, 'Color', cols(a,:));
end
xlim(ax2, [zoom2(1) zoom2(2)]);
ylim(ax2, [zoom2(3) zoom2(4)]);
set(ax2,'FontSize',10);
title(ax2,'放大区域2','FontSize',10,'FontWeight','bold');

connectRectToInset(axMain, zoom1, insetPos1);
connectRectToInset(axMain, zoom2, insetPos2);

hold(axMain,'off');
end


function zoom = resolveZoomY(zoom, iters, muAll)
x1 = max(1, round(zoom(1)));
x2 = min(numel(iters), round(zoom(2)));
if x2 <= x1, x2 = min(numel(iters), x1+50); end

y1 = zoom(3); y2 = zoom(4);
if isnan(y1) || isnan(y2)
    seg = muAll(x1:x2, :);
    mn = min(seg(:));
    mx = max(seg(:));
    pad = 0.03 * (mx - mn + eps);
    y1 = mn - pad;
    y2 = mx + pad;
end
zoom = [x1 x2 y1 y2];
end


function connectRectToInset(axMain, zoom, insetPos)
x2 = zoom(2);
y2 = zoom(4);
y1 = zoom(3);

pTR = data2norm(axMain, [x2 y2]);
pBR = data2norm(axMain, [x2 y1]);

pInsetTL = insetPos(1:2) + [0 insetPos(4)];
pInsetBL = insetPos(1:2);

annotation('line', [pTR(1) pInsetTL(1)], [pTR(2) pInsetTL(2)], 'LineWidth', 1,'Color','r');
annotation('line', [pBR(1) pInsetBL(1)], [pBR(2) pInsetBL(2)], 'LineWidth', 1,'Color','r');
end


function p = data2norm(ax, xy)
oldUnits = ax.Units;
ax.Units = 'normalized';
axPos = ax.Position;
ax.Units = oldUnits;

xl = xlim(ax); yl = ylim(ax);
px = axPos(1) + (xy(1)-xl(1))/(xl(2)-xl(1)) * axPos(3);
py = axPos(2) + (xy(2)-yl(1))/(yl(2)-yl(1)) * axPos(4);
p = [px py];
end


%% ============================================================
%  图：ECDF 对比（BestDose）
%% ============================================================
function makeECDFCompare(All, algoList, numInstances, fontName, fsAxis, fsTitle)
A = numel(algoList);
cols = lines(A);

figure('Color','w','Name','ECDF：BestDose多算法对比','Position',[140,140,980,520]);
hold on; grid on; box on;

for a = 1:A
    algoName = algoList{a}.name;
    S = All.(matlab.lang.makeValidName(algoName));
    x = S.bestDoseRuns(:);

    [f, xx] = ecdf(x);
    plot(xx, f, 'LineWidth', 2.2, 'Color', cols(a,:));
end

xlabel('BestDose（\muSv）','FontName',fontName,'FontSize',fsAxis);
ylabel('ECDF：P(BestDose \le x)','FontName',fontName,'FontSize',fsAxis);
title(sprintf('ECDF对比：%d个扰动实例 BestDose 分布', numInstances), ...
    'FontName',fontName,'FontSize',fsTitle,'FontWeight','bold');
set(gca,'FontName',fontName,'FontSize',fsAxis,'LineWidth',1.2);

leg = cell(A,1);
for a = 1:A, leg{a} = char(algoList{a}.name); end
legend(leg,'Location','southeast');

hold off;
end


%% ============================================================
%  图：速度-质量权衡散点（FirstHit vs BestDose）
%% ============================================================
function makeSpeedQualityScatter(All, algoList, numInstances, fontName, fsAxis, fsTitle)
A = numel(algoList);
cols = lines(A);

figure('Color','w','Name','速度-质量权衡：FirstHit vs BestDose','Position',[160,160,980,560]);
hold on; grid on; box on;

for a = 1:A
    algoName = algoList{a}.name;
    S = All.(matlab.lang.makeValidName(algoName));

    x = S.firstHitIter(:);
    y = S.bestDoseRuns(:);

    scatter(x, y, 40, 'o', 'MarkerFaceColor','none', ...
        'MarkerEdgeColor', cols(a,:), 'LineWidth', 1.2);
end

xlabel('FirstHitIter（迭代步，越小越快）','FontName',fontName,'FontSize',fsAxis);
ylabel('BestDose（\muSv，越小越好）','FontName',fontName,'FontSize',fsAxis);
title(sprintf('速度-质量权衡：FirstHitIter 与 BestDose（%d个实例）', numInstances), ...
    'FontName',fontName,'FontSize',fsTitle,'FontWeight','bold');
set(gca,'FontName',fontName,'FontSize',fsAxis,'LineWidth',1.2);

leg = cell(A,1);
for a = 1:A, leg{a} = char(algoList{a}.name); end
legend(leg,'Location','northeastoutside');

hold off;
end


%% ============================================================
%  图：Noise 敏感性（散点 + 线性拟合）
%% ============================================================
function makeNoiseScatterFit(All, algoList, numInstances, yField, xlab, ylab, ttl, fontName, fsAxis, fsTitle)
A = numel(algoList);
cols = lines(A);

figure('Color','w','Name',ttl,'Position',[180,180,980,560]);
hold on; grid on; box on;

for a = 1:A
    algoName = algoList{a}.name;
    S = All.(matlab.lang.makeValidName(algoName));

    x = S.instNoiseLevel(:);
    y = S.(yField)(:);

    scatter(x, y, 38, 'o', 'MarkerFaceColor','none', 'MarkerEdgeColor', cols(a,:), 'LineWidth', 1.2);

    p = polyfit(x, y, 1);
    xx = linspace(min(x), max(x), 50);
    yy = polyval(p, xx);
    plot(xx, yy, '-', 'Color', cols(a,:), 'LineWidth', 2.0);
end

xlabel(xlab,'FontName',fontName,'FontSize',fsAxis);
ylabel(ylab,'FontName',fontName,'FontSize',fsAxis);
title(ttl, 'FontName',fontName,'FontSize',fsTitle,'FontWeight','bold');
set(gca,'FontName',fontName,'FontSize',fsAxis,'LineWidth',1.2);

leg = cell(A,1);
for a = 1:A, leg{a} = char(algoList{a}.name); end
legend(leg,'Location','northeastoutside');

hold off;
end


%% ============================================================
%  图：最差情形风险（P90/P95/Max）
%% ============================================================
function makeRiskBar(All, algoList, numInstances, fontName, fsAxis, fsTitle)
A = numel(algoList);
p90 = zeros(A,1);
p95 = zeros(A,1);
mx  = zeros(A,1);

for a = 1:A
    S = All.(matlab.lang.makeValidName(algoList{a}.name));
    x = S.bestDoseRuns(:);
    p90(a) = prctile(x, 90);
    p95(a) = prctile(x, 95);
    mx(a)  = max(x);
end

figure('Color','w','Name','最差情形风险：P90/P95/Max','Position',[200,200,980,520]);
hold on; grid on; box on;

X = 1:A;
bw = 0.25;
bar(X-bw, p90, bw);
bar(X,     p95, bw);
bar(X+bw,  mx,  bw);

set(gca,'XTick',1:A,'XTickLabel',cellfun(@(s)char(s.name), algoList,'UniformOutput',false));
xlabel('算法','FontName',fontName,'FontSize',fsAxis);
ylabel('BestDose（\muSv）','FontName',fontName,'FontSize',fsAxis);
title(sprintf('最差情形风险对比：P90 / P95 / Max（%d个扰动实例）', numInstances), ...
    'FontName',fontName,'FontSize',fsTitle,'FontWeight','bold');
set(gca,'FontName',fontName,'FontSize',fsAxis,'LineWidth',1.2);

legend({'P90','P95','Max'},'Location','northeast');

hold off;
end


%% ============================================================
%  图：综合评分（归一化 + 加权）
%% ============================================================
function makeCompositeScore(All, algoList, wDose, wTime, wLen, fontName, fsAxis, fsTitle)
A = numel(algoList);

mDose = zeros(A,1);
mTime = zeros(A,1);
mLen  = zeros(A,1);

for a = 1:A
    S = All.(matlab.lang.makeValidName(algoList{a}.name));
    mDose(a) = mean(S.bestDoseRuns(:), 'omitnan');
    mTime(a) = mean(S.timeRuns(:), 'omitnan');
    mLen(a)  = mean(S.bestLenRuns(:), 'omitnan');
end

nDose = (mDose - min(mDose)) / (max(mDose) - min(mDose) + eps);
nTime = (mTime - min(mTime)) / (max(mTime) - min(mTime) + eps);
nLen  = (mLen  - min(mLen )) / (max(mLen ) - min(mLen ) + eps);

score = wDose*nDose + wTime*nTime + wLen*nLen;
[scoreSort, idx] = sort(score, 'ascend');
labels = cell(A,1);
for i = 1:A, labels{i} = char(algoList{idx(i)}.name); end

figure('Color','w','Name','综合评分（归一化加权，越小越好）','Position',[220,220,980,520]);
hold on; grid on; box on;

bar(scoreSort);
set(gca,'XTick',1:A,'XTickLabel',labels);
ylabel('综合评分（越小越好）','FontName',fontName,'FontSize',fsAxis);
title(sprintf('综合评分（min-max归一化）：wDose=%.2f, wTime=%.2f, wLen=%.2f', wDose, wTime, wLen), ...
    'FontName',fontName,'FontSize',fsTitle,'FontWeight','bold');
set(gca,'FontName',fontName,'FontSize',fsAxis,'LineWidth',1.2);

hold off;
end


%% ============================================================
%  导出：每算法摘要统计（mean/std/median/p90/p95/max）
%% ============================================================
function T = exportSummaryStats(All, algoList, outDir)
A = numel(algoList);

Algo = strings(A,1);
Dose_mean = zeros(A,1); Dose_std = zeros(A,1); Dose_med = zeros(A,1); Dose_p90 = zeros(A,1); Dose_p95 = zeros(A,1); Dose_max = zeros(A,1);
Time_mean = zeros(A,1); Time_std = zeros(A,1); Time_med = zeros(A,1);
Len_mean  = zeros(A,1); Len_std  = zeros(A,1); Len_med  = zeros(A,1);

for a = 1:A
    name = algoList{a}.name;
    S = All.(matlab.lang.makeValidName(name));
    d = S.bestDoseRuns(:);
    t = S.timeRuns(:);
    L = S.bestLenRuns(:);

    Algo(a) = name;

    Dose_mean(a) = mean(d,'omitnan'); Dose_std(a) = std(d,0,'omitnan'); Dose_med(a)=median(d,'omitnan');
    Dose_p90(a)  = prctile(d,90); Dose_p95(a)=prctile(d,95); Dose_max(a)=max(d);

    Time_mean(a) = mean(t,'omitnan'); Time_std(a)=std(t,0,'omitnan'); Time_med(a)=median(t,'omitnan');
    Len_mean(a)  = mean(L,'omitnan'); Len_std(a)=std(L,0,'omitnan'); Len_med(a)=median(L,'omitnan');
end

T = table(Algo, Dose_mean, Dose_std, Dose_med, Dose_p90, Dose_p95, Dose_max, ...
          Time_mean, Time_std, Time_med, Len_mean, Len_std, Len_med);

out = fullfile(outDir, 'summary_stats_ALL_algorithms.csv');
writetable(T, out);
fprintf('✅ 已导出统计摘要：%s\n', out);
end


%% ============================================================
%  导出：BestDose 两两 Mann-Whitney U + BH校正 + A12效应量
%% ============================================================
function T = exportPairwiseStats_BestDose(All, algoList, outDir)
A = numel(algoList);

pairs = nchoosek(1:A, 2);
P = size(pairs,1);

Algo1 = strings(P,1);
Algo2 = strings(P,1);
pval  = nan(P,1);
Ustat = nan(P,1);
A12   = nan(P,1);

for i = 1:P
    a = pairs(i,1); b = pairs(i,2);
    n1 = algoList{a}.name; n2 = algoList{b}.name;
    S1 = All.(matlab.lang.makeValidName(n1));
    S2 = All.(matlab.lang.makeValidName(n2));

    x = S1.bestDoseRuns(:);
    y = S2.bestDoseRuns(:);

    try
        [p,~,stats] = ranksum(x, y);
        pval(i) = p;
        if isfield(stats,'ranksum')
            Ustat(i) = stats.ranksum;
        end
    catch
        pval(i) = nan;
        Ustat(i) = nan;
    end

    A12(i) = vargha_delaney_A12(x, y);

    Algo1(i) = n1;
    Algo2(i) = n2;
end

p_adj = bh_fdr(pval);

T = table(Algo1, Algo2, pval, p_adj, Ustat, A12, ...
    'VariableNames', {'算法1','算法2','p值_MannWhitney','p值_BH_FDR','RankSum','A12_算法1更优概率'});

out = fullfile(outDir, 'pairwise_mannwhitney_BestDose.csv');
writetable(T, out);
fprintf('✅ 已导出两两检验：%s\n', out);
end


function A12 = vargha_delaney_A12(x, y)
x = x(:); y = y(:);
nx = numel(x); ny = numel(y);
cnt = 0; eq  = 0;
for i = 1:nx
    cnt = cnt + sum(x(i) < y);
    eq  = eq  + sum(x(i) == y);
end
A12 = (cnt + 0.5*eq) / (nx*ny + eps);
end


function p_adj = bh_fdr(p)
p_adj = nan(size(p));
idx = find(~isnan(p));
if isempty(idx), return; end
pv = p(idx);
[m, order] = sort(pv, 'ascend');
n = numel(m);
adj = zeros(n,1);
for i = 1:n
    adj(i) = m(i) * n / i;
end
for i = n-1:-1:1
    adj(i) = min(adj(i), adj(i+1));
end
adj = min(adj, 1);
tmp = zeros(n,1);
tmp(order) = adj;
p_adj(idx) = tmp;
end


%% ============================================================
%  文字结果总结：命令行输出 + txt文件
%% ============================================================
function printMergedTextSummary(All, algoList, numInstances, Tsummary, Tpairwise, ...
    outDir, doCompositeScore, wDose, wTime, wLen)

A = numel(algoList);

[~, idxDose] = sort(Tsummary.Dose_mean, 'ascend');
[~, idxTime] = sort(Tsummary.Time_mean, 'ascend');
[~, idxLen ] = sort(Tsummary.Len_mean,  'ascend');

firstHitMean = zeros(A,1);
lastImpMean  = zeros(A,1);
for a = 1:A
    name = algoList{a}.name;
    S = All.(matlab.lang.makeValidName(name));
    firstHitMean(a) = mean(S.firstHitIter(:), 'omitnan');
    lastImpMean(a)  = mean(S.lastImproveIter(:), 'omitnan');
end
[~, idxHit]  = sort(firstHitMean, 'ascend');
[~, idxLast] = sort(lastImpMean,  'ascend');

score = [];
scoreRankText = '未启用综合评分。';
if doCompositeScore
    mDose = Tsummary.Dose_mean;
    mTime = Tsummary.Time_mean;
    mLen  = Tsummary.Len_mean;

    nDose = (mDose - min(mDose)) / (max(mDose) - min(mDose) + eps);
    nTime = (mTime - min(mTime)) / (max(mTime) - min(mTime) + eps);
    nLen  = (mLen  - min(mLen )) / (max(mLen ) - min(mLen ) + eps);

    score = wDose*nDose + wTime*nTime + wLen*nLen;
    [scoreSort, idxScore] = sort(score, 'ascend');

    tmp = strings(A,1);
    for i = 1:A
        tmp(i) = sprintf('%d) %s (%.4f)', i, string(Tsummary.Algo(idxScore(i))), scoreSort(i));
    end
    scoreRankText = strjoin(cellstr(tmp), '  |  ');
else
    idxScore = [];
end

% ---------- 3) 显著性检验里找与 CPO-DBO 相关的结果 ----------
pairText = '未找到与 CPO-DBO 相关的两两检验结果。';

algo1_col = Tpairwise{:, '算法1'};
algo2_col = Tpairwise{:, '算法2'};

if isstring(algo1_col)
    algo1_col = cellstr(algo1_col);
elseif ischar(algo1_col)
    algo1_col = cellstr(algo1_col);
end

if isstring(algo2_col)
    algo2_col = cellstr(algo2_col);
elseif ischar(algo2_col)
    algo2_col = cellstr(algo2_col);
end

mask1 = strcmp(algo1_col, 'CPO-DBO');
mask2 = strcmp(algo2_col, 'CPO-DBO');
Tcmp = Tpairwise(mask1 | mask2, :);

if ~isempty(Tcmp)
    linesPair = cell(height(Tcmp),1);

    algo1_cmp = Tcmp{:, '算法1'};
    algo2_cmp = Tcmp{:, '算法2'};

    if isstring(algo1_cmp)
        algo1_cmp = cellstr(algo1_cmp);
    elseif ischar(algo1_cmp)
        algo1_cmp = cellstr(algo1_cmp);
    end

    if isstring(algo2_cmp)
        algo2_cmp = cellstr(algo2_cmp);
    elseif ischar(algo2_cmp)
        algo2_cmp = cellstr(algo2_cmp);
    end

    p_col   = Tcmp{:, 'p值_BH_FDR'};
    a12_col = Tcmp{:, 'A12_算法1更优概率'};

    for i = 1:height(Tcmp)
        a1  = algo1_cmp{i};
        a2  = algo2_cmp{i};
        p   = p_col(i);
        a12 = a12_col(i);

        if strcmp(a1, 'CPO-DBO')
            other = a2;
            probText = sprintf('A12=%.3f（CPO-DBO优于%s的概率）', a12, other);
        else
            other = a1;
            probText = sprintf('A12=%.3f（%s优于CPO-DBO的概率）', a12, other);
        end

        if isnan(p)
            sigText = 'p值缺失';
        elseif p < 0.05
            sigText = sprintf('差异显著（BH校正p=%.4g）', p);
        else
            sigText = sprintf('差异不显著（BH校正p=%.4g）', p);
        end

        linesPair{i} = sprintf('CPO-DBO vs %s: %s, %s', other, sigText, probText);
    end

    pairText = strjoin(linesPair, sprintf('\n'));
end

txt = {};
txt{end+1} = '================= 多算法汇总文字结果 =================';
txt{end+1} = sprintf('实例数：%d', numInstances);
txt{end+1} = ' ';

txt{end+1} = '一、各算法统计摘要（均值层面）';
for i = 1:height(Tsummary)
    txt{end+1} = sprintf(['%s：BestDose均值=%.4f，标准差=%.4f，P90=%.4f，P95=%.4f，Max=%.4f；' ...
        'Time均值=%.4f s；Len均值=%.4f m'], ...
        string(Tsummary.Algo(i)), ...
        Tsummary.Dose_mean(i), Tsummary.Dose_std(i), Tsummary.Dose_p90(i), ...
        Tsummary.Dose_p95(i), Tsummary.Dose_max(i), ...
        Tsummary.Time_mean(i), Tsummary.Len_mean(i));
end
txt{end+1} = ' ';

txt{end+1} = '二、按不同指标排序';
txt{end+1} = sprintf('1) 最优总剂量均值（越小越好）：%s', joinRankText(Tsummary.Algo, Tsummary.Dose_mean, idxDose));
txt{end+1} = sprintf('2) 运行时间均值（越小越好）：%s', joinRankText(Tsummary.Algo, Tsummary.Time_mean, idxTime));
txt{end+1} = sprintf('3) 路径总距离均值（越小越好）：%s', joinRankText(Tsummary.Algo, Tsummary.Len_mean, idxLen));
txt{end+1} = sprintf('4) FirstHitIter均值（越小越快）：%s', joinRankText(Tsummary.Algo, firstHitMean, idxHit));
txt{end+1} = sprintf('5) LastImproveIter均值（越小越早稳定）：%s', joinRankText(Tsummary.Algo, lastImpMean, idxLast));
txt{end+1} = ' ';

txt{end+1} = '三、综合评分';
txt{end+1} = scoreRankText;
txt{end+1} = ' ';

txt{end+1} = '四、与 CPO-DBO 相关的显著性检验';
txt{end+1} = char(pairText);
txt{end+1} = ' ';

bestDoseAlgo = string(Tsummary.Algo(idxDose(1)));
bestTimeAlgo = string(Tsummary.Algo(idxTime(1)));
bestLenAlgo  = string(Tsummary.Algo(idxLen(1)));

if doCompositeScore
    bestCompAlgo = string(Tsummary.Algo(idxScore(1)));
    txt{end+1} = '五、自动总结';
    txt{end+1} = sprintf(['从20个扰动实例的统计结果来看，%s在最优总剂量均值方面表现最佳，' ...
        '%s在运行时间方面占优，%s在路径总距离方面表现最佳。综合考虑剂量、时间和距离后，' ...
        '%s的综合评分最低，可视为总体表现最优的算法。'], ...
        bestDoseAlgo, bestTimeAlgo, bestLenAlgo, bestCompAlgo);
else
    txt{end+1} = '五、自动总结';
    txt{end+1} = sprintf(['从20个扰动实例的统计结果来看，%s在最优总剂量均值方面表现最佳，' ...
        '%s在运行时间方面占优，%s在路径总距离方面表现最佳。'], ...
        bestDoseAlgo, bestTimeAlgo, bestLenAlgo);
end

txt{end+1} = '=======================================================';

fprintf('\n');
for i = 1:numel(txt)
    fprintf('%s\n', txt{i});
end

outTxt = fullfile(outDir, 'text_summary_ALL_algorithms.txt');
fid = fopen(outTxt, 'w', 'n', 'UTF-8');
if fid ~= -1
    for i = 1:numel(txt)
        fprintf(fid, '%s\n', txt{i});
    end
    fclose(fid);
    fprintf('✅ 已导出文字总结：%s\n', outTxt);
else
    warning('无法写入文字总结文件：%s', outTxt);
end

end


function s = joinRankText(algoNames, vals, idx)
parts = strings(numel(idx),1);
for i = 1:numel(idx)
    k = idx(i);
    parts(i) = sprintf('%d) %s (%.4f)', i, string(algoNames(k)), vals(k));
end
s = strjoin(cellstr(parts), '  |  ');
end


%% ===================== utils =====================
function name = getFileName(fullpath)
[~, name, ext] = fileparts(fullpath);
name = [name ext];
end