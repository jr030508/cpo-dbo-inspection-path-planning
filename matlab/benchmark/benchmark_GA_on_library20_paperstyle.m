function benchmark_GA_on_library20_paperstyle()
% benchmark_GA_on_library20_paperstyle.m
% ============================================================
%  GA 接入母本（Step2+3） + 论文风统计与绘图（与前面风格一致）
%
%  输入：
%   - instance_lib20/instance_library20.mat
%       必含：doseMatrices(NxNxK), nodeData, xCoords, yCoords, numInstances
%       可选：instNoiseLevel, instStrengths, instX, instY
%
%  输出（默认 ./bench_GA_20_paper/ ）：
%   - ranked_results_GA_20instances.csv
%   - results_GA_20instances.mat
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
    instX = S.instX; instY = S.instY;
else
    instX = nan(numInstances,7); instY = nan(numInstances,7);
end

%% ===================== 1) 输出目录 =====================
outDir = fullfile(pwd,'bench_GA_20_paper');
if ~exist(outDir,'dir'), mkdir(outDir); end

fprintf('\n================= GA：母本20实例公平跑批 + 论文风统计绘图 =================\n');
fprintf('母本：%s\n', libMat);
fprintf('输出：%s\n', outDir);

%% ===================== 2) GA 参数（与你原版一致） =====================
opts = struct();
opts.popSize   = 400;
opts.maxGen    = 5000;
opts.pc        = 0.85;
opts.pm        = 0.15;
opts.eliteRatio = 0.05;
opts.eliteNum  = max(1, round(opts.popSize * opts.eliteRatio));

opts.tournamentSize = 3;

opts.baseSeed = 7;

%% ===================== 3) 跑批（20实例） =====================
K = numInstances;
maxIter = opts.maxGen;

bestDoseRuns = zeros(K,1);
bestLenRuns  = zeros(K,1);
timeRuns     = zeros(K,1);
firstHitIter = zeros(K,1);
lastImproveIter = zeros(K,1);
sceneFP      = zeros(K,1);

bestPathCell = cell(K,1);
convAll = nan(maxIter, K);

AllResults = struct();
AllResults.algoName = "GA";
AllResults.opts = opts;
AllResults.libraryMat = libMat;

for k = 1:K
    rng(opts.baseSeed + k, 'twister');
    doseMatrix = doseMatrices(:,:,k);

    tStart = tic;
    result = run_GA_oneInstance(doseMatrix, numNodes, opts);
    tCost = toc(tStart);

    bestDoseRuns(k) = result.bestDose;
    timeRuns(k)     = tCost;
    bestPathCell{k} = result.bestPath;

    bestLenRuns(k)  = pathLengthMeters(result.bestPath, xCoords, yCoords);
    convAll(:,k)    = result.convergence(:);

    [firstHitIter(k), lastImproveIter(k)] = computeHitImprove(result.convergence, result.bestDose);

    if hasXY && hasStrength && hasNoise
        sources_k  = [instX(k,:).', instY(k,:).'];
        strength_k = instStrengths(k,:);
        noise_k    = instNoiseLevel(k);
        sceneFP(k) = makeSceneFP(sources_k, strength_k, noise_k);
    else
        sceneFP(k) = makeDoseFP(doseMatrix);
    end

    writematrix(result.bestPath(:).', fullfile(outDir, sprintf('bestPath_inst%02d.csv', k)));

    fprintf('Inst %02d/%02d | BestDose=%.4f | BestLen=%.1f | Hit=%d | LastImp=%d | Time=%.2fs | Noise=%.3f | FP=%.6f\n', ...
        k, K, bestDoseRuns(k), bestLenRuns(k), firstHitIter(k), lastImproveIter(k), timeRuns(k), instNoiseLevel(k), sceneFP(k));

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

%% ===================== 4) 排名表 + CSV =====================
[bestDoseSorted, idxSort] = sort(bestDoseRuns, 'ascend');
instSorted  = idxSort(:);

bestLenSorted = bestLenRuns(instSorted);
hitSorted     = firstHitIter(instSorted);
lastImpSorted = lastImproveIter(instSorted);
timeSorted    = timeRuns(instSorted);
noiseSorted   = instNoiseLevel(instSorted);
fpSorted      = sceneFP(instSorted);

fprintf('\n================= GA：按最优剂量排序（从低到高）=================\n');
fprintf('Rank | Inst | BestDose | BestLen | FirstHit | LastImprove | Time(s) | Noise | sceneFP\n');
for r = 1:K
    fprintf('%4d | %4d | %8.4f | %7.1f | %8d | %11d | %7.3f | %.3f | %.6f\n', ...
        r, instSorted(r), bestDoseSorted(r), bestLenSorted(r), ...
        hitSorted(r), lastImpSorted(r), timeSorted(r), noiseSorted(r), fpSorted(r));
end

T = table((1:K)', instSorted(:), bestDoseSorted(:), bestLenSorted(:), ...
          hitSorted(:), lastImpSorted(:), timeSorted(:), noiseSorted(:), fpSorted(:), ...
          'VariableNames', {'Rank','Instance','BestDose','BestLen','FirstHitIter','LastImproveIter','TimeSec','NoiseLevel','SceneFP'});

rankCsv = fullfile(outDir, 'ranked_results_GA_20instances.csv');
writetable(T, rankCsv);
fprintf('\n已导出：%s\n', rankCsv);

%% ===================== 5) 论文风图（同前） =====================
iters = (1:maxIter)';

TopK = 5;
Ktop = min(TopK, K);
topIdx = instSorted(1:Ktop);
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
    '最优总剂量（\muSv）', sprintf('GA：20个扰动实例最优总剂量（箱线图+散点，Top-%d高亮）', Ktop), ...
    fontName, fsAxis, fsTitle, lwBox, lwMedian, jit, mkSize, mkLW, grayEdge, grayFace, topColor);

plotBoxScatter(bestLenRuns, isTop, Ktop, ...
    '最优路径总距离（m）', sprintf('GA：20个扰动实例最优路径总距离（箱线图+散点，Top-%d高亮）', Ktop), ...
    fontName, fsAxis, fsTitle, lwBox, lwMedian, jit, mkSize, mkLW, grayEdge, grayFace, topColor);

plotBoxScatter(timeRuns, isTop, Ktop, ...
    '运行时间（s）', sprintf('GA：20个扰动实例运行时间（箱线图+散点，Top-%d高亮）', Ktop), ...
    fontName, fsAxis, fsTitle, lwBox, lwMedian, jit, mkSize, mkLW, grayEdge, grayFace, topColor);

convMean = mean(convAll, 2, 'omitnan');
convStd  = std(convAll, 0, 2, 'omitnan');

figMean = figure('Color','w','Name','收敛曲线：均值±标准差（GA）','Position',[160,160,1060,480]);
hold on; grid on; box on;
xShade = [iters; flipud(iters)];
yShade = [convMean-convStd; flipud(convMean+convStd)];
hFill = fill(xShade, yShade, shadeColor, 'EdgeColor','none');
hFill.FaceAlpha = shadeAlpha;
hMean = plot(iters, convMean, '-', 'LineWidth', lwMean);
hMean.Color = meanColor;

xlabel('迭代次数','FontName',fontName,'FontSize',fsAxis);
ylabel('最优剂量（\muSv）','FontName',fontName,'FontSize',fsAxis);
title(sprintf('GA 收敛曲线（均值±1标准差）：%d个扰动实例（maxGen=%d）', K, maxIter), ...
    'FontName',fontName,'FontSize',fsTitle,'FontWeight','bold');
set(gca,'FontName',fontName,'FontSize',fsAxis,'LineWidth',1.2);
legend([hFill, hMean], {'±1σ 区间','均值曲线'}, 'Location','northeast');
hold off;

figTop = figure('Color','w','Name',sprintf('Top-%d 收敛曲线叠加（GA，含局部放大）',Ktop),'Position',[180,180,1120,520]);
hold on; grid on; box on;
cols = lines(Ktop);
for kk = 1:Ktop
    rInst = topIdx(kk);
    plot(iters, convAll(:,rInst), 'LineWidth', lwCurve, 'Color', cols(kk,:));
end
xlabel('迭代次数','FontName',fontName,'FontSize',fsAxis);
ylabel('最优剂量（\muSv）','FontName',fontName,'FontSize',fsAxis);
title(sprintf('GA：Top-%d 最优实例收敛曲线叠加（按最优剂量排序）', Ktop), ...
    'FontName',fontName,'FontSize',fsTitle,'FontWeight','bold');
set(gca,'FontName',fontName,'FontSize',fsAxis,'LineWidth',1.2);

leg = cell(Ktop,1);
for kk = 1:Ktop
    rInst = topIdx(kk);
    leg{kk} = sprintf('Top%02d：实例=%d，最优=%.3f，hit=%d', kk, rInst, bestDoseRuns(rInst), firstHitIter(rInst));
end
legend(leg,'Location','northeastoutside');

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

%% ===================== 6) 保存结果 =====================
save(fullfile(outDir,'results_GA_20instances.mat'), ...
    'bestDoseRuns','bestLenRuns','timeRuns','convAll','maxIter','numInstances', ...
    'firstHitIter','lastImproveIter','idxSort','topIdx','sceneFP', ...
    'instNoiseLevel','instStrengths','instX','instY','opts','AllResults');

fprintf('\n✅ 完成：GA（20实例）统计与论文风绘图已生成。\n');
fprintf('✅ 已保存：%s\n', fullfile(outDir,'results_GA_20instances.mat'));

saveas(figMean, fullfile(outDir, sprintf('GA_mean_std_conv_%dinst.png',K)));
saveas(figTop,  fullfile(outDir, sprintf('GA_top%d_conv_inset.png',Ktop)));

end % ===== end main =====


%% ============================================================
%                    GA: one instance (permutation-safe)
%% ============================================================
function result = run_GA_oneInstance(doseMatrix, numNodes, opts)

popSize = opts.popSize;
maxGen  = opts.maxGen;
pc      = opts.pc;
pm      = opts.pm;
eliteNum = opts.eliteNum;
tSize   = opts.tournamentSize;

% init pop
population = cell(popSize,1);
fitness = zeros(popSize,1);
for i = 1:popSize
    population{i} = [1, randperm(numNodes-1)+1, 1];
    fitness(i) = calculateFitness(population{i}, doseMatrix);
end

[bestFitness, bestIdx] = min(fitness);
bestPath = population{bestIdx};

convergence = zeros(maxGen,1);

for gen = 1:maxGen
    % selection
    newPop = selectionTournament(population, fitness, popSize, tSize);

    % crossover (OX)
    newPop = crossoverOX(newPop, pc);

    % mutation (swap + occasional inversion)
    newPop = mutationSwap(newPop, pm);

    % eval new
    newFit = zeros(popSize,1);
    for i = 1:popSize
        newPop{i} = validatePath_det(newPop{i}, numNodes);
        newFit(i) = calculateFitness(newPop{i}, doseMatrix);
    end

    % elite from old
    [sortedFitOld, idxOld] = sort(fitness, 'ascend');
    elitePaths = population(idxOld(1:eliteNum));
    eliteFit   = sortedFitOld(1:eliteNum);

    % merge & keep best popSize
    allPop = [newPop; elitePaths];
    allFit = [newFit; eliteFit];

    [sortedAllFit, idxAll] = sort(allFit, 'ascend');
    population = allPop(idxAll(1:popSize));
    fitness    = sortedAllFit(1:popSize);

    if fitness(1) < bestFitness
        bestFitness = fitness(1);
        bestPath = population{1};
    end

    convergence(gen) = bestFitness;

    if mod(gen, 200) == 0
        fprintf('  GA gen %d/%d | BestDose=%.6f\n', gen, maxGen, bestFitness);
    end
end

bestIter = find(convergence <= bestFitness + 1e-12, 1, 'first');
if isempty(bestIter), bestIter = find(convergence==min(convergence), 1, 'first'); end

result.bestDose = bestFitness;
result.bestPath = bestPath;
result.bestIter = bestIter;
result.convergence = convergence;

end


%% ============================================================
%                    GA operators (permutation-safe)
%% ============================================================
function newPop = selectionTournament(pop, fit, popSize, tSize)
newPop = cell(popSize,1);
N = numel(pop);
for i = 1:popSize
    cand = randperm(N, tSize);
    [~, idx] = min(fit(cand));
    newPop{i} = pop{cand(idx)};
end
end

function newPop = crossoverOX(pop, pc)
newPop = pop;
for i = 1:2:numel(pop)-1
    if rand() < pc
        p1 = pop{i};   p2 = pop{i+1};
        c1 = oxCrossoverFull(p1, p2);
        c2 = oxCrossoverFull(p2, p1);
        newPop{i}   = c1;
        newPop{i+1} = c2;
    end
end
end

function child = oxCrossoverFull(parentA, parentB)
a = parentA(2:end-1);
b = parentB(2:end-1);
c = oxPerm(a, b);
child = [1, c, 1];
end

function c = oxPerm(a, b)
n = numel(a);
c = zeros(1,n);
i = randi([1, n-1]);
j = randi([i+1, n]);
c(i:j) = a(i:j);
fillPos = [1:i-1, j+1:n];
bSeq = b(~ismember(b, c(i:j)));
c(fillPos) = bSeq;
end

function newPop = mutationSwap(pop, pm)
newPop = pop;
for i = 1:numel(pop)
    if rand() < pm
        p = pop{i};
        mid = p(2:end-1);
        n = numel(mid);
        idx = randperm(n,2);
        tmp = mid(idx(1)); mid(idx(1)) = mid(idx(2)); mid(idx(2)) = tmp;

        % 10%概率加一次 inversion，增强跳出局部（仍保持排列）
        if rand() < 0.10
            a = randi([1, n-1]); b = randi([a+1, n]);
            mid(a:b) = fliplr(mid(a:b));
        end
        newPop{i} = [1, mid, 1];
    end
end
end

function totalDose = calculateFitness(path, doseMatrix)
totalDose = 0;
for i = 1:(length(path)-1)
    totalDose = totalDose + doseMatrix(path(i), path(i+1));
end
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