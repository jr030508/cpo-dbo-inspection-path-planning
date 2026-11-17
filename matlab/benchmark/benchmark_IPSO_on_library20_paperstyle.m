function benchmark_IPSO_on_library20_paperstyle()
% benchmark_IPSO_on_library20_paperstyle.m
% ============================================================
%  IPSO 接入母本（Step2+3） + 论文风统计与绘图（与你给的CPO风格一致）
%
%  输入：
%   - instance_lib20/instance_library20.mat
%       必含：doseMatrices(NxNxK), nodeData, xCoords, yCoords, numInstances
%       可选：instNoiseLevel, instStrengths, instX, instY
%
%  输出（默认 ./bench_IPSO_20_paper/ ）：
%   - ranked_results_IPSO_20instances.csv
%   - results_IPSO_20instances.mat
%   - bestPath_instXX.csv
%   - 论文风图：箱线图+散点 / 均值±std / TopK叠加+放大窗
% ============================================================

clc; close all;

%% ===================== 0) 读母本 =====================
libMat = fullfile(pwd,'instance_lib20','instance_library20.mat');
assert(exist(libMat,'file')==2, '找不到母本：%s', libMat);
S = load(libMat);

req = {'doseMatrices','nodeData','xCoords','yCoords','numInstances'};
for i=1:numel(req)
    assert(isfield(S,req{i}), '母本缺字段：%s', req{i});
end

doseMatrices = S.doseMatrices;
nodeData     = S.nodeData;
xCoords      = S.xCoords;
yCoords      = S.yCoords;
numInstances = S.numInstances;

numNodes = size(nodeData,1);
assert(ndims(doseMatrices)==3 && size(doseMatrices,3)>=numInstances, 'doseMatrices 维度应为 NxNxK');

% 可选字段（用于 Noise / SceneFP）
hasNoise    = isfield(S,'instNoiseLevel');
hasStrength = isfield(S,'instStrengths');
hasXY       = isfield(S,'instX') && isfield(S,'instY');

if hasNoise,    instNoiseLevel = S.instNoiseLevel(:); else, instNoiseLevel = nan(numInstances,1); end
if hasStrength, instStrengths  = S.instStrengths;     else, instStrengths  = nan(numInstances,7); end
if hasXY
    instX = S.instX; instY = S.instY;  % 20x7
else
    instX = nan(numInstances,7); instY = nan(numInstances,7);
end

%% ===================== 1) 输出目录 =====================
outDir = fullfile(pwd,'bench_IPSO_20_paper');
if ~exist(outDir,'dir'), mkdir(outDir); end

fprintf('\n================= IPSO：母本20实例公平跑批 + 论文风统计绘图 =================\n');
fprintf('母本：%s\n', libMat);
fprintf('输出：%s\n', outDir);

%% ===================== 2) IPSO 参数 =====================
opts = struct();
opts.M            = 600;
opts.maxChaosIter = 200;
opts.maxPSOIter   = 5000;
opts.maxIterTotal = opts.maxChaosIter + opts.maxPSOIter;

% permutation-safe 变异概率（一次PSO更新后触发一次）
opts.p_swap      = 0.60;
opts.p_insert    = 0.15;
opts.p_inversion = 0.15;
opts.p_scramble  = 0.10;

% OX交叉中“局部用gBest片段”的概率（保留你原0.1思想）
opts.p_use_gBest_segment = 0.10;

% 每实例固定种子（公平+可复现）
opts.baseSeed = 7;

%% ===================== 3) 跑批（20实例） =====================
K = numInstances;
maxIter = opts.maxIterTotal;

bestDoseRuns = zeros(K,1);
bestLenRuns  = zeros(K,1);
timeRuns     = zeros(K,1);
firstHitIter = zeros(K,1);
lastImproveIter = zeros(K,1);
sceneFP      = zeros(K,1);

bestPathCell = cell(K,1);
convAll = nan(maxIter, K);  % 每列一个实例的收敛曲线（长度 maxIter）

for k = 1:K
    rng(opts.baseSeed + k, 'twister');
    doseMatrix = doseMatrices(:,:,k);

    tStart = tic;
    [gBestPath, gBestFitness, conv] = IPSO_oneInstance(doseMatrix, numNodes, opts);
    tCost = toc(tStart);

    % 路径长度（几何距离，仅用于报告，不参与优化）
    bestLen = pathLengthMeters(gBestPath, xCoords, yCoords);

    bestDoseRuns(k) = gBestFitness;
    bestLenRuns(k)  = bestLen;
    timeRuns(k)     = tCost;
    bestPathCell{k} = gBestPath;

    % 统一 convAll
    L = min(numel(conv), maxIter);
    convAll(1:L, k) = conv(1:L);

    % first hit / last improve（按该实例的最终最优值）
    [firstHitIter(k), lastImproveIter(k)] = computeHitImprove(conv, gBestFitness);

    % sceneFP（如果母本有源项/强度/噪声就做真实指纹；否则用doseMatrix做指纹）
    if hasXY && hasStrength && hasNoise
        sources_k = [instX(k,:).', instY(k,:).']; % 7x2
        strength_k = instStrengths(k,:);
        noise_k = instNoiseLevel(k);
        sceneFP(k) = makeSceneFP(sources_k, strength_k, noise_k);
    else
        sceneFP(k) = makeDoseFP(doseMatrix);
    end

    % 保存路径
    writematrix(gBestPath(:).', fullfile(outDir, sprintf('bestPath_inst%02d.csv', k)));

    fprintf('Inst %02d/%02d | BestDose=%.4f | BestLen=%.1f | Hit=%d | LastImp=%d | Time=%.2fs | Noise=%.3f | FP=%.6f\n', ...
        k, K, bestDoseRuns(k), bestLenRuns(k), firstHitIter(k), lastImproveIter(k), timeRuns(k), instNoiseLevel(k), sceneFP(k));
end

%% ===================== 4) 按最优剂量排序 + 打印排名表 + 导出CSV =====================
[bestDoseSorted, idxSort] = sort(bestDoseRuns, 'ascend');
instSorted  = idxSort(:);

bestLenSorted = bestLenRuns(instSorted);
hitSorted     = firstHitIter(instSorted);
lastImpSorted = lastImproveIter(instSorted);
timeSorted    = timeRuns(instSorted);
noiseSorted   = instNoiseLevel(instSorted);
fpSorted      = sceneFP(instSorted);

fprintf('\n================= IPSO：按最优暴露代价排序（从低到高）=================\n');
fprintf('Rank | Inst | BestDose | BestLen | FirstHit | LastImprove | Time(s) | Noise | sceneFP\n');
for r = 1:K
    fprintf('%4d | %4d | %8.4f | %7.1f | %8d | %11d | %7.3f | %.3f | %.6f\n', ...
        r, instSorted(r), bestDoseSorted(r), bestLenSorted(r), ...
        hitSorted(r), lastImpSorted(r), timeSorted(r), noiseSorted(r), fpSorted(r));
end

T = table((1:K)', instSorted(:), bestDoseSorted(:), bestLenSorted(:), ...
          hitSorted(:), lastImpSorted(:), timeSorted(:), noiseSorted(:), fpSorted(:), ...
          'VariableNames', {'Rank','Instance','BestDose','BestLen','FirstHitIter','LastImproveIter','TimeSec','NoiseLevel','SceneFP'});

rankCsv = fullfile(outDir, 'ranked_results_IPSO_20instances.csv');
writetable(T, rankCsv);
fprintf('\n已导出：%s\n', rankCsv);

%% ===================== 5) 论文风图：箱线图+散点 / 均值±std / TopK叠加+放大窗 =====================
iters = (1:maxIter)';

TopK = 5;
Ktop = min(TopK, K);
topIdx = instSorted(1:Ktop);         % 这里 topIdx 是“实例编号”
isTop = false(K,1);
isTop(topIdx) = true;                % 注意：实例编号=1..K，所以可直接索引

theme = "orange"; % "orange" / "purpleblue"
switch theme
    case "orange"
        shadeColor = [0.95, 0.45, 0.10];
        meanColor  = [0.05, 0.05, 0.05];
        topColor   = [0.90, 0.25, 0.05];
    case "purpleblue"
        shadeColor = [0.35, 0.45, 0.95];
        meanColor  = [0.05, 0.05, 0.05];
        topColor   = [0.20, 0.35, 0.95];
end

shadeAlpha = 0.30;
lwMean     = 2.4;
lwCurve    = 2.2;

fontName = 'SimSun';
fsAxis   = 12;
fsTitle  = 14;

lwBox    = 1.8;
lwMedian = 2.4;
jit      = 0.10;
mkSize   = 40;
mkLW     = 1.4;
grayEdge = [0.25 0.25 0.25];
grayFace = [0.20 0.20 0.20];

plotBoxScatter(bestDoseRuns, isTop, Ktop, ...
    '最优总暴露代价（a.u.）', sprintf('IPSO：20个扰动实例最优总暴露代价（箱线图+散点，Top-%d高亮）', Ktop), ...
    fontName, fsAxis, fsTitle, lwBox, lwMedian, jit, mkSize, mkLW, grayEdge, grayFace, topColor);

plotBoxScatter(bestLenRuns, isTop, Ktop, ...
    '最优路径总距离（m）', sprintf('IPSO：20个扰动实例最优路径总距离（箱线图+散点，Top-%d高亮）', Ktop), ...
    fontName, fsAxis, fsTitle, lwBox, lwMedian, jit, mkSize, mkLW, grayEdge, grayFace, topColor);

plotBoxScatter(timeRuns, isTop, Ktop, ...
    '运行时间（s）', sprintf('IPSO：20个扰动实例运行时间（箱线图+散点，Top-%d高亮）', Ktop), ...
    fontName, fsAxis, fsTitle, lwBox, lwMedian, jit, mkSize, mkLW, grayEdge, grayFace, topColor);

% 均值±1σ（按行：迭代维度）
convMean = mean(convAll, 2, 'omitnan');
convStd  = std(convAll, 0, 2, 'omitnan');

figMean = figure('Color','w','Name','收敛曲线：均值±标准差（IPSO）','Position',[160,160,1060,480]);
hold on; grid on; box on;

xShade = [iters; flipud(iters)];
yShade = [convMean-convStd; flipud(convMean+convStd)];
hFill = fill(xShade, yShade, shadeColor, 'EdgeColor','none');
hFill.FaceAlpha = shadeAlpha;

hMean = plot(iters, convMean, '-', 'LineWidth', lwMean);
hMean.Color = meanColor;

xline(opts.maxChaosIter, 'r--', 'LineWidth', 1.2, 'Label', '混沌->PSO');

xlabel('迭代次数','FontName',fontName,'FontSize',fsAxis);
ylabel('最优暴露代价（a.u.）','FontName',fontName,'FontSize',fsAxis);
title(sprintf('IPSO 收敛曲线（均值±1标准差）：%d个扰动实例（Chaos=%d, PSO=%d）', K, opts.maxChaosIter, opts.maxPSOIter), ...
    'FontName',fontName,'FontSize',fsTitle,'FontWeight','bold');
set(gca,'FontName',fontName,'FontSize',fsAxis,'LineWidth',1.2);
legend([hFill, hMean], {'±1σ 区间','均值曲线'}, 'Location','northeast');
hold off;

% Top-K 收敛曲线叠加 + inset 放大
figTop = figure('Color','w','Name',sprintf('Top-%d 收敛曲线叠加（IPSO，含局部放大）',Ktop),'Position',[180,180,1120,520]);
hold on; grid on; box on;

cols = lines(Ktop);
for kk = 1:Ktop
    rInst = topIdx(kk);              % 实例编号
    plot(iters, convAll(:, rInst), 'LineWidth', lwCurve, 'Color', cols(kk,:));
end
xline(opts.maxChaosIter, 'r--', 'LineWidth', 1.2, 'Label', '混沌->PSO');

xlabel('迭代次数','FontName',fontName,'FontSize',fsAxis);
ylabel('最优暴露代价（a.u.）','FontName',fontName,'FontSize',fsAxis);
title(sprintf('IPSO：Top-%d 最优实例收敛曲线叠加（按最优暴露代价排序）', Ktop), ...
    'FontName',fontName,'FontSize',fsTitle,'FontWeight','bold');
set(gca,'FontName',fontName,'FontSize',fsAxis,'LineWidth',1.2);

leg = cell(Ktop,1);
for kk = 1:Ktop
    rInst = topIdx(kk);
    leg{kk} = sprintf('Top%02d：实例=%d，最优=%.3f，hit=%d', kk, rInst, bestDoseRuns(rInst), firstHitIter(rInst));
end
legend(leg,'Location','northeastoutside');

% inset zoom（默认放大前600次迭代，可按你习惯改）
x1 = 0;  x2 = 600;
topEnd = min(bestDoseRuns(topIdx));
y2 = topEnd * 1.05;
y1 = topEnd * 0.95;
insetPos = [0.20 0.38 0.22 0.30];

axMain = gca;
rectangle(axMain, 'Position', [x1, y1, x2-x1, y2-y1], ...
    'EdgeColor', [0 0 0], 'LineWidth', 1.5, 'LineStyle', '-');

axInset = axes('Position', insetPos);
box(axInset, 'on'); grid(axInset, 'on'); hold(axInset, 'on');

linesMain = findobj(axMain, 'Type', 'line');
linesMain = flipud(linesMain); % 保持顺序
for kk = 1:numel(linesMain)
    copyobj(linesMain(kk), axInset);
end
xlim(axInset, [x1 x2]);
ylim(axInset, [y1 y2]);
set(axInset, 'FontSize', 10);
title(axInset, '局部放大', 'FontSize', 10);

p1 = data2norm(axMain, [x2 y2]);
p2 = insetPos(1:2) + [0 insetPos(4)];
annotation('line', [p1(1) p2(1)], [p1(2) p2(2)], 'LineWidth', 1);

p1b = data2norm(axMain, [x2 y1]);
p2b = insetPos(1:2);
annotation('line', [p1b(1) p2b(1)], [p1b(2) p2b(2)], 'LineWidth', 1);

hold off;

%% ===================== 6) 保存结果（与你CPO风格一致） =====================
save(fullfile(outDir,'results_IPSO_20instances.mat'), ...
    'bestDoseRuns','bestLenRuns','timeRuns','convAll','maxIter','numInstances', ...
    'firstHitIter','lastImproveIter','idxSort','topIdx','sceneFP', ...
    'instNoiseLevel','instStrengths','instX','instY','opts');

fprintf('\n✅ 完成：IPSO（20实例）统计与论文风绘图已生成。\n');
fprintf('✅ 已保存：%s\n', fullfile(outDir,'results_IPSO_20instances.mat'));

% 同时把关键图保存为png（便于论文插图）
saveas(figMean, fullfile(outDir, sprintf('IPSO_mean_std_conv_%dinst.png',K)));
saveas(figTop,  fullfile(outDir, sprintf('IPSO_top%d_conv_inset.png',Ktop)));

end % ===== end main =====


%% ============================================================
%                   IPSO core (Permutation-safe)
%% ============================================================
function [gBestPath, gBestFitness, convergence] = IPSO_oneInstance(doseMatrix, numNodes, opts)

M            = opts.M;
maxChaosIter = opts.maxChaosIter;
maxPSOIter   = opts.maxPSOIter;

dim = numNodes - 1;                    % 对应节点2..N
chaosPop = initChaosPopulation(M, dim);

population = cell(M,1);
fitness    = zeros(M,1);

for i = 1:M
    path = chaosToPath_permSafe(chaosPop(i,:), numNodes);
    population{i} = path;
    fitness(i) = calculateFitness(path, doseMatrix);
end

pBestPath    = population;
pBestFitness = fitness;
[gBestFitness, gBestIdx] = min(fitness);
gBestPath = population{gBestIdx};

convergence = zeros(maxChaosIter + maxPSOIter, 1);

% ---- chaos stage
for iter = 1:maxChaosIter
    chaosPop = updateChaos(chaosPop);

    for i = 1:M
        newPath = chaosToPath_permSafe(chaosPop(i,:), numNodes);
        newFit  = calculateFitness(newPath, doseMatrix);
        if newFit < pBestFitness(i)
            pBestFitness(i) = newFit;
            pBestPath{i}    = newPath;
        end
    end

    [minFit, minIdx] = min(pBestFitness);
    if minFit < gBestFitness
        gBestFitness = minFit;
        gBestPath    = pBestPath{minIdx};
    end

    convergence(iter) = gBestFitness;
end

% ---- PSO stage (F1/F2/F3)
for iter = 1:maxPSOIter
    for i = 1:M
        cur = population{i};

        % F1: OX(cur, pBest)
        c1 = oxCrossover(cur, pBestPath{i});

        % F2: OX(c1, gBest)
        c2 = oxCrossover(c1, gBestPath);

        % (可选) 10%概率直接替换一段为gBest片段（保留你原思想）
        if rand() < opts.p_use_gBest_segment
            c2 = injectSegmentFromGBest(c2, gBestPath);
        end

        % F3: mutation
        newPath = mutate_permSafe(c2, numNodes, opts);

        % deterministic validate
        newPath = validatePath_det(newPath, numNodes);

        newFit = calculateFitness(newPath, doseMatrix);

        if newFit < pBestFitness(i)
            pBestFitness(i) = newFit;
            pBestPath{i}    = newPath;
        end

        population{i} = newPath;
        fitness(i)    = newFit;
    end

    [minFit, minIdx] = min(pBestFitness);
    if minFit < gBestFitness
        gBestFitness = minFit;
        gBestPath    = pBestPath{minIdx};
    end

    convergence(maxChaosIter + iter) = gBestFitness;
end

end

function chaosPop = initChaosPopulation(M, dim)
chaosPop = zeros(M, dim);
badVals = [0, 0.25, 0.5, 0.75, 1];
for i = 1:M
    temp = rand(1, dim);
    while any(ismembertol(temp, badVals, 1e-12))
        temp = rand(1, dim);
    end
    chaosPop(i,:) = temp;
end
end

function newChaos = updateChaos(oldChaos)
newChaos = 4 .* oldChaos .* (1 - oldChaos);
end

function path = chaosToPath_permSafe(chaosVector, numNodes)
[~, idx] = sort(chaosVector, 'ascend');
perm = idx + 1; % 2..N
path = [1, perm, 1];
end

function childPath = oxCrossover(parentA, parentB)
a = parentA(2:end-1);
b = parentB(2:end-1);
c = oxPerm(a, b);
childPath = [1, c, 1];
end

function child = oxPerm(a, b)
n = numel(a);
child = zeros(1,n);
i = randi([1, n-1]);
j = randi([i+1, n]);
child(i:j) = a(i:j);
fillPos = [1:i-1, j+1:n];
bSeq = b(~ismember(b, child(i:j)));
child(fillPos) = bSeq;
end

function p2 = mutate_permSafe(path, numNodes, opts)
p = path(2:end-1);
if rand() < opts.p_swap,      p = swapPerm(p); end
if rand() < opts.p_insert,    p = insertPerm(p); end
if rand() < opts.p_inversion, p = inversionPerm(p); end
if rand() < opts.p_scramble,  p = scramblePerm(p); end
p2 = [1, p, 1];
p2 = validatePath_det(p2, numNodes);
end

function p = swapPerm(p)
n = numel(p);
idx = randperm(n,2);
tmp = p(idx(1)); p(idx(1)) = p(idx(2)); p(idx(2)) = tmp;
end

function p = insertPerm(p)
n = numel(p);
i = randi([1,n]); j = randi([1,n]);
while j == i, j = randi([1,n]); end
val = p(i);
p(i) = [];
p = [p(1:j-1), val, p(j:end)];
end

function p = inversionPerm(p)
n = numel(p);
i = randi([1, n-1]); j = randi([i+1, n]);
p(i:j) = fliplr(p(i:j));
end

function p = scramblePerm(p)
n = numel(p);
i = randi([1, n-1]); j = randi([i+1, n]);
seg = p(i:j);
p(i:j) = seg(randperm(numel(seg)));
end

function validPath = validatePath_det(path, numNodes)
path(1) = 1; path(end) = 1;
mid = path(2:end-1);
mid = mid(mid>=2 & mid<=numNodes);
[~, ia] = unique(mid, 'stable');
mid = mid(ia);
allNodes = 2:numNodes;
missing = setdiff(allNodes, mid, 'stable');
mid = [mid, missing];
mid = mid(1:(numNodes-1));
validPath = [1, mid, 1];
end

function totalDose = calculateFitness(path, doseMatrix)
totalDose = 0;
for i = 1:(length(path)-1)
    totalDose = totalDose + doseMatrix(path(i), path(i+1));
end
end

function path = injectSegmentFromGBest(path, gBest)
% 在不破坏排列的情况下，强制把某段改成 gBest 的对应段（然后修复）
p = path(2:end-1);
g = gBest(2:end-1);
n = numel(p);

i = randi([1, n-1]);
j = randi([i+1, n]);

p(i:j) = g(i:j);

% 修复重复/缺失：按出现顺序保留，再补齐缺失
[~, ia] = unique(p, 'stable');
p = p(sort(ia));
missing = setdiff(2:(n+1), p, 'stable'); % 2..numNodes
p = [p, missing];
p = p(1:n);

path = [1, p, 1];
end

%% ============================================================
%                 Metrics helpers (Len/Hit/FP)
%% ============================================================
function L = pathLengthMeters(path, x, y)
L = 0;
for i = 1:(numel(path)-1)
    a = path(i); b = path(i+1);
    L = L + hypot(x(a)-x(b), y(a)-y(b));
end
end

function [firstHit, lastImp] = computeHitImprove(conv, finalBest)
epsTol = 1e-12;
hitIdx = find(conv <= finalBest + epsTol, 1, 'first');
if isempty(hitIdx), hitIdx = numel(conv); end
firstHit = hitIdx;

% last improve：最后一次严格下降的位置
d = diff(conv);
impIdx = find(d < -1e-12);
if isempty(impIdx)
    lastImp = 1;
else
    lastImp = impIdx(end) + 1;
end
end

function fp = makeSceneFP(sources, strength, noiseLevel)
fp = noiseLevel + 1e-6*sum(strength(:)) + 1e-3*sum(sources(:));
end

function fp = makeDoseFP(doseMatrix)
% 轻量指纹：用于确认同一实例（不依赖源项字段）
% 取若干统计量混合，保证稳定
v = doseMatrix(:);
v(v==0) = [];
p1 = prctile(v, 50);
p2 = prctile(v, 90);
fp = 1e-6*sum(v) + 1e-3*p2 + 1e-4*p1;
end

%% ============================================================
%                 Plot helpers (与你给的风格一致)
%% ============================================================
function plotBoxScatter(vec, isTop, K, ylab, ttl, fontName, fsAxis, fsTitle, lwBox, lwMedian, jit, mkSize, mkLW, grayEdge, grayFace, topColor)
vec = vec(:);
n = numel(vec);

figure('Color','w','Name',ttl,'Position',[100,100,980,420]);
hold on; grid on; box on;

boxplot(vec, 'Labels',{''}, 'Symbol','', 'Whisker',1.5);
set(findobj(gca,'Tag','Box'),     'LineWidth', lwBox);
set(findobj(gca,'Tag','Whisker'), 'LineWidth', lwBox);
set(findobj(gca,'Tag','Cap'),     'LineWidth', lwBox);
set(findobj(gca,'Tag','Median'),  'LineWidth', lwMedian);

hBox = findobj(gca,'Tag','Box');
for kk = 1:numel(hBox)
    xd = get(hBox(kk),'XData'); yd = get(hBox(kk),'YData');
    patch(xd, yd, grayFace, 'FaceAlpha',0.10, 'EdgeColor', grayEdge, 'LineWidth', lwBox);
end

x0 = 1;
xj = x0 + (rand(n,1)-0.5)*2*jit;

scatter(xj(~isTop), vec(~isTop), mkSize, 'o', ...
    'MarkerFaceColor','none', 'MarkerEdgeColor', grayEdge, 'LineWidth', mkLW);

scatter(xj(isTop), vec(isTop), mkSize*1.2, 'p', ...
    'MarkerFaceColor','none', 'MarkerEdgeColor', topColor, 'LineWidth', mkLW+0.3);

ylabel(ylab,'FontName',fontName,'FontSize',fsAxis);
title(ttl, 'FontName',fontName,'FontSize',fsTitle,'FontWeight','bold');
set(gca,'FontName',fontName,'FontSize',fsAxis,'LineWidth',1.2);
xlim([0.6 1.4]); set(gca,'XTick',[]);

legend({'','普通实例',sprintf('Top-%d实例',K)},'Location','northeast');
hold off;
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

function name = getFileName(fullpath)
[~, name, ext] = fileparts(fullpath);
name = [name ext];
end