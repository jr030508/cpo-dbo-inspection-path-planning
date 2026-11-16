function benchmark_CPO_DBO_on_library20_paperstyle()
% benchmark_CPO_DBO_on_library20_paperstyle.m
% ============================================================
%  CPO-DBO 接入母本（Step2+3） + 论文风统计与绘图（与你给的CPO风格一致）
%
%  输入：
%   - instance_lib20/instance_library20.mat
%       必含：doseMatrices(NxNxK), nodeData, xCoords, yCoords, numInstances
%       可选：instNoiseLevel, instStrengths, instX, instY
%
%  输出（默认 ./bench_CPO_DBO_20_paper/ ）：
%   - ranked_results_CPO_DBO_20instances.csv
%   - results_CPO_DBO_20instances.mat
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
outDir = fullfile(pwd,'bench_CPO_DBO_20_paper');
if ~exist(outDir,'dir'), mkdir(outDir); end

fprintf('\n================= CPO-DBO：母本20实例公平跑批 + 论文风统计绘图 =================\n');
fprintf('母本：%s\n', libMat);
fprintf('输出：%s\n', outDir);

%% ===================== 2) CPO-DBO 参数（与你原版一致） =====================
opts = struct();

opts.M        = 400;
opts.maxIter  = 5000;
opts.eliteKeep = max(2, round(0.02*opts.M));

opts.Tscale   = 30000;

opts.popBall  = round(opts.M * 0.4);
opts.popBreed = round(opts.M * 0.2);
opts.popBaby  = opts.M - opts.popBall - opts.popBreed;

opts.defenseProb0 = 0.70;
opts.attackProb0  = 0.25;
opts.escapeProb0  = 0.05;

opts.twoOpt_defense = 1;
opts.twoOpt_attack  = 1;
opts.twoOpt_escape  = 0;

opts.mut_defense = [0.60, 0.15, 0.15];
opts.mut_attack  = [0.35, 0.20, 0.30];
opts.mut_escape  = [0.25, 0.25, 0.45];

opts.baseSeed = 7;

%% ===================== 3) 跑批（20实例） =====================
K = numInstances;
maxIter = opts.maxIter;

bestDoseRuns = zeros(K,1);
bestLenRuns  = zeros(K,1);
timeRuns     = zeros(K,1);
firstHitIter = zeros(K,1);
lastImproveIter = zeros(K,1);
sceneFP      = zeros(K,1);

bestPathCell = cell(K,1);
convAll = nan(maxIter, K);  % 每列一个实例的收敛曲线

AllResults = struct();
AllResults.algoName = "CPO-DBO";
AllResults.opts = opts;
AllResults.libraryMat = libMat;

for k = 1:K
    rng(opts.baseSeed + k, 'twister');
    doseMatrix = doseMatrices(:,:,k);

    tStart = tic;
    result = run_CPO_DBO_oneInstance(doseMatrix, numNodes, opts);
    tCost = toc(tStart);

    bestDoseRuns(k) = result.bestDose;
    timeRuns(k)     = tCost;
    bestPathCell{k} = result.bestPath;

    % 几何距离（仅用于报告）
    bestLenRuns(k)  = pathLengthMeters(result.bestPath, xCoords, yCoords);

    % 收敛曲线
    convAll(:,k) = result.convergence(:);

    % first hit / last improve
    [firstHitIter(k), lastImproveIter(k)] = computeHitImprove(result.convergence, result.bestDose);

    % sceneFP
    if hasXY && hasStrength && hasNoise
        sources_k = [instX(k,:).', instY(k,:).']; % 7x2
        strength_k = instStrengths(k,:);
        noise_k = instNoiseLevel(k);
        sceneFP(k) = makeSceneFP(sources_k, strength_k, noise_k);
    else
        sceneFP(k) = makeDoseFP(doseMatrix);
    end

    % 保存路径
    writematrix(result.bestPath(:).', fullfile(outDir, sprintf('bestPath_inst%02d.csv', k)));

    fprintf('Inst %02d/%02d | BestDose=%.4f | BestLen=%.1f | Hit=%d | LastImp=%d | Time=%.2fs | Noise=%.3f | FP=%.6f\n', ...
        k, K, bestDoseRuns(k), bestLenRuns(k), firstHitIter(k), lastImproveIter(k), timeRuns(k), instNoiseLevel(k), sceneFP(k));

    % 保存到结构体（保留原有信息）
    AllResults.instances(k).instanceId  = k;
    AllResults.instances(k).bestDose    = result.bestDose;
    AllResults.instances(k).bestIter    = result.bestIter;
    AllResults.instances(k).runtimeS    = tCost;
    AllResults.instances(k).bestPath    = result.bestPath;
    AllResults.instances(k).convergence = result.convergence;

    if hasNoise,    AllResults.instances(k).noiseLevel = instNoiseLevel(k); end
    if hasStrength, AllResults.instances(k).strengths  = instStrengths(k,:); end
    if hasXY,       AllResults.instances(k).sourcesXY  = [instX(k,:).', instY(k,:).']; end
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

fprintf('\n================= CPO-DBO：按最优剂量排序（从低到高）=================\n');
fprintf('Rank | Inst | BestDose | BestLen | FirstHit | LastImprove | Time(s) | Noise | sceneFP\n');
for r = 1:K
    fprintf('%4d | %4d | %8.4f | %7.1f | %8d | %11d | %7.3f | %.3f | %.6f\n', ...
        r, instSorted(r), bestDoseSorted(r), bestLenSorted(r), ...
        hitSorted(r), lastImpSorted(r), timeSorted(r), noiseSorted(r), fpSorted(r));
end

T = table((1:K)', instSorted(:), bestDoseSorted(:), bestLenSorted(:), ...
          hitSorted(:), lastImpSorted(:), timeSorted(:), noiseSorted(:), fpSorted(:), ...
          'VariableNames', {'Rank','Instance','BestDose','BestLen','FirstHitIter','LastImproveIter','TimeSec','NoiseLevel','SceneFP'});

rankCsv = fullfile(outDir, 'ranked_results_CPO_DBO_20instances.csv');
writetable(T, rankCsv);
fprintf('\n已导出：%s\n', rankCsv);

%% ===================== 5) 论文风图：箱线图+散点 / 均值±std / TopK叠加+放大窗 =====================
iters = (1:maxIter)';

TopK = 5;
Ktop = min(TopK, K);
topIdx = instSorted(1:Ktop);     % topIdx 是“实例编号”
isTop = false(K,1);
isTop(topIdx) = true;

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
    '最优总剂量（\muSv）', sprintf('CPO-DBO：20个扰动实例最优总剂量（箱线图+散点，Top-%d高亮）', Ktop), ...
    fontName, fsAxis, fsTitle, lwBox, lwMedian, jit, mkSize, mkLW, grayEdge, grayFace, topColor);

plotBoxScatter(bestLenRuns, isTop, Ktop, ...
    '最优路径总距离（m）', sprintf('CPO-DBO：20个扰动实例最优路径总距离（箱线图+散点，Top-%d高亮）', Ktop), ...
    fontName, fsAxis, fsTitle, lwBox, lwMedian, jit, mkSize, mkLW, grayEdge, grayFace, topColor);

plotBoxScatter(timeRuns, isTop, Ktop, ...
    '运行时间（s）', sprintf('CPO-DBO：20个扰动实例运行时间（箱线图+散点，Top-%d高亮）', Ktop), ...
    fontName, fsAxis, fsTitle, lwBox, lwMedian, jit, mkSize, mkLW, grayEdge, grayFace, topColor);

% 均值±1σ
convMean = mean(convAll, 2, 'omitnan');
convStd  = std(convAll, 0, 2, 'omitnan');

figMean = figure('Color','w','Name','收敛曲线：均值±标准差（CPO-DBO）','Position',[160,160,1060,480]);
hold on; grid on; box on;

xShade = [iters; flipud(iters)];
yShade = [convMean-convStd; flipud(convMean+convStd)];
hFill = fill(xShade, yShade, shadeColor, 'EdgeColor','none');
hFill.FaceAlpha = shadeAlpha;

hMean = plot(iters, convMean, '-', 'LineWidth', lwMean);
hMean.Color = meanColor;

xlabel('迭代次数','FontName',fontName,'FontSize',fsAxis);
ylabel('最优剂量（\muSv）','FontName',fontName,'FontSize',fsAxis);
title(sprintf('CPO-DBO 收敛曲线（均值±1标准差）：%d个扰动实例（maxIter=%d）', K, maxIter), ...
    'FontName',fontName,'FontSize',fsTitle,'FontWeight','bold');
set(gca,'FontName',fontName,'FontSize',fsAxis,'LineWidth',1.2);
legend([hFill, hMean], {'±1σ 区间','均值曲线'}, 'Location','northeast');
hold off;

% Top-K 收敛叠加 + inset
figTop = figure('Color','w','Name',sprintf('Top-%d 收敛曲线叠加（CPO-DBO，含局部放大）',Ktop),'Position',[180,180,1120,520]);
hold on; grid on; box on;
cols = lines(Ktop);
for kk = 1:Ktop
    rInst = topIdx(kk); % 实例编号
    plot(iters, convAll(:,rInst), 'LineWidth', lwCurve, 'Color', cols(kk,:));
end
xlabel('迭代次数','FontName',fontName,'FontSize',fsAxis);
ylabel('最优剂量（\muSv）','FontName',fontName,'FontSize',fsAxis);
title(sprintf('CPO-DBO：Top-%d 最优实例收敛曲线叠加（按最优剂量排序）', Ktop), ...
    'FontName',fontName,'FontSize',fsTitle,'FontWeight','bold');
set(gca,'FontName',fontName,'FontSize',fsAxis,'LineWidth',1.2);

leg = cell(Ktop,1);
for kk = 1:Ktop
    rInst = topIdx(kk);
    leg{kk} = sprintf('Top%02d：实例=%d，最优=%.3f，hit=%d', kk, rInst, bestDoseRuns(rInst), firstHitIter(rInst));
end
legend(leg,'Location','northeastoutside');

% inset zoom（默认放大前600次迭代）
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
linesMain = flipud(linesMain);
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

%% ===================== 6) 保存结果（论文可直接用） =====================
save(fullfile(outDir,'results_CPO_DBO_20instances.mat'), ...
    'bestDoseRuns','bestLenRuns','timeRuns','convAll','maxIter','numInstances', ...
    'firstHitIter','lastImproveIter','idxSort','topIdx','sceneFP', ...
    'instNoiseLevel','instStrengths','instX','instY','opts','AllResults');

fprintf('\n✅ 完成：CPO-DBO（20实例）统计与论文风绘图已生成。\n');
fprintf('✅ 已保存：%s\n', fullfile(outDir,'results_CPO_DBO_20instances.mat'));

% 同时保存关键图
saveas(figMean, fullfile(outDir, sprintf('CPO_DBO_mean_std_conv_%dinst.png',K)));
saveas(figTop,  fullfile(outDir, sprintf('CPO_DBO_top%d_conv_inset.png',Ktop)));

end % ===== end main =====


%% ========================================================================
%                          CPO-DBO 单实例主函数（保持原逻辑）
% ========================================================================
function result = run_CPO_DBO_oneInstance(doseMatrix, numNodes, opts)

M       = opts.M;
maxIter = opts.maxIter;
eliteKeep = opts.eliteKeep;

popBall  = opts.popBall;
popBreed = opts.popBreed;

defenseProb0 = opts.defenseProb0;
attackProb0  = opts.attackProb0;
escapeProb0  = opts.escapeProb0;

Tscale = opts.Tscale;

twoOptIters_defense = opts.twoOpt_defense;
twoOptIters_attack  = opts.twoOpt_attack;
twoOptIters_escape  = opts.twoOpt_escape;

mut_defense = opts.mut_defense;
mut_attack  = opts.mut_attack;
mut_escape  = opts.mut_escape;

population = cell(M,1);
fitness    = zeros(M,1);

for i = 1:M
    perm = randperm(numNodes-1) + 1;
    population{i} = [1, perm, 1];
    fitness(i) = calculateFitness(population{i}, doseMatrix);
end

pBestPath    = population;
pBestFitness = fitness;

[gBestFitness, gBestIdx] = min(fitness);
gBestPath = population{gBestIdx};

convergence = zeros(maxIter, 1);

for iter = 1:maxIter
    t = min(1, iter / Tscale);

    defenseProb = defenseProb0 * (1 - 0.6*t);
    attackProb  = attackProb0  * (1 - 0.2*t);
    escapeProb  = min(0.60, escapeProb0 + 0.55*t);

    sProb = defenseProb + attackProb + escapeProb;
    defenseProb = defenseProb / sProb;
    attackProb  = attackProb  / sProb;
    escapeProb  = escapeProb  / sProb;

    [~, idxSort] = sort(fitness, 'ascend');
    eliteIdx = idxSort(1:eliteKeep);

    newPopulation = population;
    newFitness    = fitness;

    for e = 1:eliteKeep
        newPopulation{e} = population{eliteIdx(e)};
        newFitness(e)    = fitness(eliteIdx(e));
    end

    for i = eliteKeep+1:M
        currentPath = population{i};
        currentFit  = fitness(i);

        if i <= popBall
            groupBias = 0.10;
        elseif i <= popBall + popBreed
            groupBias = -0.05;
        else
            groupBias = -0.10;
        end

        r = rand();
        pD = min(0.95, max(0.05, defenseProb + groupBias));
        pA = min(0.95, max(0.05, attackProb  - groupBias/2));
        pS = pD + pA + escapeProb;
        pD = pD / pS;
        pA = pA / pS;

        if r < pD
            newPath = defenseBehavior_perm(currentPath, gBestPath, numNodes);
            newPath = mutate_perm(newPath, numNodes, mut_defense(1), mut_defense(2), mut_defense(3));
            newPath = twoOptLocal(newPath, doseMatrix, twoOptIters_defense);

        elseif r < pD + pA
            mateIdx  = randi(M);
            matePath = pBestPath{mateIdx};

            newPath = oxCrossover(currentPath, matePath);
            newPath = mutate_perm(newPath, numNodes, mut_attack(1), mut_attack(2), mut_attack(3));
            newPath = twoOptLocal(newPath, doseMatrix, twoOptIters_attack);

        else
            newPath = escapeBehavior_perm(currentPath, numNodes);
            newPath = mutate_perm(newPath, numNodes, mut_escape(1), mut_escape(2), mut_escape(3));
            newPath = twoOptLocal(newPath, doseMatrix, twoOptIters_escape);
        end

        newPath = validatePath_det(newPath, numNodes);
        newFit  = calculateFitness(newPath, doseMatrix);

        if newFit < currentFit
            newPopulation{i} = newPath;
            newFitness(i)    = newFit;

            if newFit < pBestFitness(i)
                pBestPath{i} = newPath;
                pBestFitness(i) = newFit;
            end
        else
            newPopulation{i} = currentPath;
            newFitness(i)    = currentFit;
        end
    end

    population = newPopulation;
    fitness    = newFitness;

    [minFitness, minIdx] = min(fitness);
    if minFitness < gBestFitness
        gBestFitness = minFitness;
        gBestPath    = population{minIdx};
    end

    convergence(iter) = gBestFitness;

    if mod(iter, 200) == 0
        fprintf('  Iter %d/%d | BestDose=%.6f | pD=%.2f pA=%.2f pE=%.2f\n', ...
            iter, maxIter, gBestFitness, defenseProb, attackProb, escapeProb);
    end
end

% bestIter：首次达到最终最优
idx = find(convergence <= gBestFitness + 1e-12, 1, 'first');
if isempty(idx), idx = find(convergence==min(convergence), 1, 'first'); end

result.bestDose    = gBestFitness;
result.bestPath    = gBestPath;
result.bestIter    = idx;
result.convergence = convergence;

end


%% ========================================================================
%                     CPO-DBO permutation-safe ops
% ========================================================================
function newPath = defenseBehavior_perm(path, bestPath, numNodes)
p = path(2:end-1);
b = bestPath(2:end-1);

if rand() < 0.80
    child = oxPerm(p, b);
else
    child = p;
end

if rand() < 0.50
    child = swapPerm(child);
end

newPath = [1, child, 1];
newPath = validatePath_det(newPath, numNodes);
end

function newPath = escapeBehavior_perm(path, numNodes)
p = path(2:end-1);

if rand() < 0.80
    p = scramblePerm(p);
end
if rand() < 0.60
    p = inversionPerm(p);
end

newPath = [1, p, 1];
newPath = validatePath_det(newPath, numNodes);
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

function path2 = mutate_perm(path, numNodes, pSwap, pInsert, pInv)
p = path(2:end-1);
if rand() < pSwap,   p = swapPerm(p); end
if rand() < pInsert, p = insertPerm(p); end
if rand() < pInv,    p = inversionPerm(p); end
path2 = [1, p, 1];
path2 = validatePath_det(path2, numNodes);
end

function p = swapPerm(p)
n = numel(p);
idx = randperm(n,2);
tmp = p(idx(1));
p(idx(1)) = p(idx(2));
p(idx(2)) = tmp;
end

function p = insertPerm(p)
n = numel(p);
i = randi([1,n]);
j = randi([1,n]);
while j == i
    j = randi([1,n]);
end
val = p(i);
p(i) = [];
p = [p(1:j-1), val, p(j:end)];
end

function p = inversionPerm(p)
n = numel(p);
i = randi([1, n-1]);
j = randi([i+1, n]);
p(i:j) = fliplr(p(i:j));
end

function p = scramblePerm(p)
n = numel(p);
i = randi([1, n-1]);
j = randi([i+1, n]);
seg = p(i:j);
p(i:j) = seg(randperm(numel(seg)));
end

function path2 = twoOptLocal(path, doseMatrix, iters)
if iters <= 0
    path2 = path;
    return;
end
p = path;
n = numel(p);
for k = 1:iters
    i = randi([2, n-3]);
    j = randi([i+1, n-2]);
    if j <= i, continue; end
    old = doseMatrix(p(i-1), p(i)) + doseMatrix(p(j), p(j+1));
    neu = doseMatrix(p(i-1), p(j)) + doseMatrix(p(i), p(j+1));
    if neu < old
        p(i:j) = fliplr(p(i:j));
    end
end
path2 = p;
end

function totalDose = calculateFitness(path, doseMatrix)
totalDose = 0;
for i = 1:(length(path)-1)
    totalDose = totalDose + doseMatrix(path(i), path(i+1));
end
end

function validPath = validatePath_det(path, numNodes)
path(1) = 1;
path(end) = 1;
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

%% ========================================================================
%                 Metrics helpers (Len/Hit/FP)
% ========================================================================
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
v = doseMatrix(:);
v(v==0) = [];
p1 = prctile(v, 50);
p2 = prctile(v, 90);
fp = 1e-6*sum(v) + 1e-3*p2 + 1e-4*p1;
end

%% ========================================================================
%                 Plot helpers (与你给的风格一致)
% ========================================================================
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