%% ============================================================
%  Case3 (51 nodes) 核场景巡检路径规划（无障碍-精细模型）
%  全部边均采用 LOS 直线 + 沿线剂量积分
%
%  目的：
%    - 作为第三部分“含障碍混合边模型”的无障碍对照基线
%    - 与含障碍版本保持统一剂量计算口径
%
%  建模方式：
%    - 所有节点对之间默认可直连（无障碍）
%    - 边长度 = 欧氏距离
%    - 边剂量 = 沿直线对辐射场采样积分（trapz）
%    - 边总代价 = Dose + epsLen * Length
%
%  优化算法：
%    - 新 CPO-DBO（Permutation-safe）
%    - OX交叉 + swap/insert/inversion(+scramble用于escape)
%    - 动态概率 + 分组偏置 + 精英保留 + 轻量2-opt
%
%  输出文件：
%    - cost_matrix_LOS_case3.csv
%    - dose_matrix_LOS_case3.csv
%    - pathlen_matrix_LOS_case3.csv
%    - edgeType_matrix_LOS_case3.csv   (全0，表示全部为LOS直线边)
%
%  出图：
%    1) 收敛曲线
%    2) 场景图（节点+点源，不含路径）
%    3) 最优路径图（全部LOS直线）
%    4) 剂量分解条形图
%% ============================================================

clc; clear; close all;
rng(6);

%% ===================== 0) USER PARAMETERS ===================
% ---- Primary objective: DOSE ----
% epsLen is a tiny tie-break for length
epsLen = 1e-5;

% ---- Straight-line dose integration step (meters) ----
dsLine = 1.0;

% ---- CPO-DBO parameters ----
M = 400;
maxIter = 5000;

eliteKeep = max(2, round(0.052*M));   % 精英保留
twoOptIters_defense = 1;
twoOptIters_attack  = 1;
twoOptIters_escape  = 0;

% DBO group split
popBall  = round(M * 0.4);
popBreed = round(M * 0.2);
popBaby  = M - popBall - popBreed; %#ok<NASGU>

% dynamic behavior probabilities
defenseProb0 = 0.70;
attackProb0  = 0.25;
escapeProb0  = 0.05;

%% ===================== 1) Nodes & sources ===================
nodeData = [
    1, 120, 0;
    2, 131, 18;
    3, 154, 11;
    4, 164, 24;
    5, 174, 35;
    6, 193, 48;
    7, 187, 67;
    8, 165, 73;
    9, 110, 16;
    10, 73, 12;
    11, 55, 29;
    12, 72, 46;
    13, 33, 52;
    14, 13, 64;
    15, 29, 81;
    16, 62, 79;
    17, 175, 91;
    18, 198, 91;
    19, 189, 99;
    20, 173, 106;
    21, 53, 92;
    22, 69, 104;
    23, 98, 107;
    24, 61, 119;
    25, 87, 127;
    26, 186, 118;
    27, 172, 151;
    28, 160, 139;
    29, 156, 150;
    30, 143, 123;
    31, 131, 134;
    32, 43, 135;
    33, 12, 143;
    34, 52, 152;
    35, 38, 161;
    36, 22, 172;
    37, 40, 192;
    38, 123, 143;
    39, 135, 166;
    40, 140, 172;
    41, 140, 183;
    42, 109, 173;
    43, 96, 193;
    44, 87, 179;
    45, 73, 166;
    46, 64, 188;
    47, 83, 25;
    48, 42, 52;
    49, 18, 88;
    50, 23, 156;
    51, 117, 181;
];
xCoords = nodeData(:,2);
yCoords = nodeData(:,3);
numNodes = size(nodeData,1);

sources = [
    35, 35;
    30, 110;
    88, 70;
    88, 142;
    138, 107;
    160, 50;
    168, 180;
];
sourceStrength = [8000, 12000, 10000, 15000, 9000, 13000, 11000];

%% ===================== 2) Radiation field for visualization ===
xMin = 0; xMax = 200;
yMin = 0; yMax = 200;

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

%% ===================== 3) Build LOS pairwise matrices =========
doseMat  = inf(numNodes, numNodes);
lenMat   = inf(numNodes, numNodes);
costMat  = inf(numNodes, numNodes);
edgeType = zeros(numNodes, numNodes); % 0=LOS straight only

fprintf('\n[STEP] 构建无障碍-精细模型边矩阵：全部采用LOS直线积分...\n');
totalPairs = numNodes*(numNodes-1)/2;
cnt = 0;

for i = 1:numNodes
    for j = i+1:numNodes
        cnt = cnt + 1;

        p1 = [xCoords(i), yCoords(i)];
        p2 = [xCoords(j), yCoords(j)];

        L = hypot(p2(1)-p1(1), p2(2)-p1(2));
        D = lineDoseIntegral(p1, p2, sources, sourceStrength, dsLine);

        lenMat(i,j)  = L; lenMat(j,i)  = L;
        doseMat(i,j) = D; doseMat(j,i) = D;

        costMat(i,j) = D + epsLen * L;
        costMat(j,i) = costMat(i,j);

        edgeType(i,j) = 0;
        edgeType(j,i) = 0;

        if mod(cnt,100)==0
            fprintf('  LOS edges: %d/%d pairs done\n', cnt, totalPairs);
        end
    end
end

for i = 1:numNodes
    lenMat(i,i) = 0;
    doseMat(i,i) = 0;
    costMat(i,i) = 0;
    edgeType(i,i) = 0;
end

writematrix(lenMat,  'pathlen_matrix_LOS_case3.csv');
writematrix(doseMat, 'dose_matrix_LOS_case3.csv');
writematrix(costMat, 'cost_matrix_LOS_case3.csv');
writematrix(edgeType,'edgeType_matrix_LOS_case3.csv');

fprintf('\nSaved:\n- pathlen_matrix_LOS_case3.csv\n- dose_matrix_LOS_case3.csv\n- cost_matrix_LOS_case3.csv\n- edgeType_matrix_LOS_case3.csv\n');
fprintf('Settings: epsLen=%.2e, dsLine=%.2f\n', epsLen, dsLine);

%% ===================== 4) 新 CPO-DBO on LOS cost matrix ======
population = cell(M,1);
fitness    = zeros(M,1);

for i = 1:M
    perm = randperm(numNodes-1) + 1;
    population{i} = [1, perm, 1];
    fitness(i) = tourCost(population{i}, costMat);
end

pBestPath = population;
pBestFitness = fitness;
[gBestFitness, gBestIdx] = min(fitness);
gBestPath = population{gBestIdx};

convergence = zeros(maxIter,1);

fprintf('\n[STEP] Start 新CPO-DBO on LOS cost matrix (Dose + epsLen*Len)...\n');

for iter = 1:maxIter
    t = iter / maxIter;

    defenseProb = defenseProb0 * (1 - 0.6*t);
    attackProb  = attackProb0  * (1 - 0.2*t);
    escapeProb  = min(0.8, escapeProb0 + 0.55*t);

    sProb = defenseProb + attackProb + escapeProb;
    defenseProb = defenseProb / sProb;
    attackProb  = attackProb  / sProb;
    escapeProb  = escapeProb  / sProb;

    % 精英保留
    [~, idxSort] = sort(fitness, 'ascend');
    newPopulation = population;
    newFitness = fitness;

    if eliteKeep > 0
        eliteIdx = idxSort(1:eliteKeep);
        for e = 1:eliteKeep
            newPopulation{e} = population{eliteIdx(e)};
            newFitness(e) = fitness(eliteIdx(e));
        end
    end

    startIdx = eliteKeep + 1;

    for i = startIdx:M
        currentPath = population{i};
        currentFit  = fitness(i);

        % 分组偏置
        if i <= popBall
            groupBias = 0.10;
        elseif i <= popBall + popBreed
            groupBias = -0.05;
        else
            groupBias = -0.10;
        end

        % 选择行为
        r = rand();
        pD = min(0.95, max(0.05, defenseProb + groupBias));
        pA = min(0.95, max(0.05, attackProb  - groupBias/2));
        pS = pD + pA + escapeProb;
        pD = pD / pS;
        pA = pA / pS;
        pE = escapeProb / pS; %#ok<NASGU>

        if r < pD
            newPath = defenseBehavior_perm(currentPath, gBestPath, numNodes);
            newPath = mutate_perm(newPath, numNodes, 0.60, 0.15, 0.15);
            newPath = twoOptLocal(newPath, costMat, twoOptIters_defense);

        elseif r < pD + pA
            mateIdx = randi(M);
            matePath = pBestPath{mateIdx};
            newPath = oxCrossover(currentPath, matePath);
            newPath = mutate_perm(newPath, numNodes, 0.35, 0.20, 0.30);
            newPath = twoOptLocal(newPath, costMat, twoOptIters_attack);

        else
            newPath = escapeBehavior_perm(currentPath, numNodes);
            newPath = mutate_perm(newPath, numNodes, 0.25, 0.25, 0.45);
            newPath = twoOptLocal(newPath, costMat, twoOptIters_escape);
        end

        newPath = validatePath_det(newPath, numNodes);
        newFit  = tourCost(newPath, costMat);

        if newFit < currentFit
            newPopulation{i} = newPath;
            newFitness(i) = newFit;

            if newFit < pBestFitness(i)
                pBestPath{i} = newPath;
                pBestFitness(i) = newFit;
            end
        else
            newPopulation{i} = currentPath;
            newFitness(i) = currentFit;
        end
    end

    population = newPopulation;
    fitness = newFitness;

    [minFit, minIdx] = min(fitness);
    if minFit < gBestFitness
        gBestFitness = minFit;
        gBestPath = population{minIdx};
    end

    convergence(iter) = gBestFitness;

    if mod(iter,50)==0
        fprintf('Iter %d/%d | BestCost=%.6f | (pD=%.2f,pA=%.2f,pE=%.2f)\n', ...
            iter, maxIter, gBestFitness, defenseProb, attackProb, escapeProb);
    end
end

%% ===================== 5) Output summary =====================
totalDose = tourCost(gBestPath, doseMat);
totalLen  = tourCost(gBestPath, lenMat);
nLOSEdges = numel(gBestPath) - 1;

fprintf('\n===== Case3 无障碍-精细模型（全LOS直线积分）+ 新CPO-DBO =====\n');
fprintf('epsLen=%.2e (tie-break), dsLine=%.2f\n', epsLen, dsLine);
fprintf('Best tour COST (Dose + epsLen*Len): %.6f\n', gBestFitness);
fprintf('Total Dose (main): %.6f\n', totalDose);
fprintf('Total Length (ref): %.6f\n', totalLen);
fprintf('LOS edges in best tour: %d / %d\n', nLOSEdges, numNodes);
fprintf('Best tour: %s\n', strjoin(cellstr(num2str(gBestPath')), '->'));

%% ===================== 6) Visualization (中文) ================
% 6.1 收敛曲线
figure('Name','收敛曲线（无障碍-精细模型+新CPO-DBO）','Position',[100,100,900,420]);
plot(convergence,'LineWidth',2); grid on;
xlabel('迭代次数');
ylabel('最优巡检代价（剂量 + \epsilon \cdot 距离）');
title('新CPO-DBO 在无障碍-精细模型边权矩阵上的收敛曲线');

% 6.2 场景图（不含路径）
figure('Name','场景示意图（无障碍-精细模型，不含路径）','Position',[100,100,1000,820]);
hold on; axis equal; grid on;
xlim([xMin xMax]); ylim([yMin yMax]);
xlabel('X坐标'); ylabel('Y坐标');
title('核场景地图：节点与点源（无障碍-精细模型）');

contourf(XX, YY, radGrid, 30, 'LineStyle','none');
colormap(jet); cb = colorbar; cb.Label.String = '辐射水平（相对值）';

scatter(xCoords, yCoords, 70, 'k', 'filled', 'MarkerEdgeColor','w');
for i = 1:numNodes
    text(xCoords(i)+2, yCoords(i)+1, sprintf('%d',i), ...
        'Color','w','FontWeight','bold', 'BackgroundColor',[0 0 0 0.35], ...
        'HorizontalAlignment','center');
end

scatter(sources(:,1), sources(:,2), 160, 'r', 'pentagram', 'filled');
for s = 1:size(sources,1)
    text(sources(s,1), sources(s,2)-5, sprintf('源%d',s), ...
        'Color','r','FontWeight','bold', 'HorizontalAlignment','center');
end

hLegNode = plot(nan,nan,'o','MarkerFaceColor','k','MarkerEdgeColor','w','LineStyle','none');
hLegSrc  = plot(nan,nan,'p','MarkerFaceColor','r','MarkerEdgeColor','r','LineStyle','none');
legend([hLegNode, hLegSrc], {'巡检节点','点源'}, 'Location','northeastoutside');
hold off;

% 6.3 最优路径图
figure('Name','最优巡检路径（无障碍-精细模型+新CPO-DBO）','Position',[100,100,1150,880]);
hold on; axis equal; grid on;
xlim([xMin xMax]); ylim([yMin yMax]);
xlabel('X坐标'); ylabel('Y坐标');
title(sprintf('最优巡检回路（总剂量=%.3f，总距离=%.1f，全部LOS，\\epsilon=%.1e）', ...
    totalDose, totalLen, epsLen));

contourf(XX, YY, radGrid, 30, 'LineStyle','none');
colormap(jet); cb = colorbar; cb.Label.String = '辐射水平（相对值）';

scatter(xCoords, yCoords, 70, 'k', 'filled', 'MarkerEdgeColor','w');
for i = 1:numNodes
    text(xCoords(i)+2, yCoords(i)+1, sprintf('%d',i), ...
        'Color','w','FontWeight','bold', 'BackgroundColor',[0 0 0 0.35], ...
        'HorizontalAlignment','center');
end

scatter(sources(:,1), sources(:,2), 160, 'r', 'pentagram', 'filled');
for s = 1:size(sources,1)
    text(sources(s,1), sources(s,2)-5, sprintf('源%d',s), ...
        'Color','r','FontWeight','bold', 'HorizontalAlignment','center');
end

cols = lines(length(gBestPath)-1);
for k = 1:(length(gBestPath)-1)
    a = gBestPath(k);
    b = gBestPath(k+1);

    p1 = [xCoords(a), yCoords(a)];
    p2 = [xCoords(b), yCoords(b)];

    plot([p1(1) p2(1)], [p1(2) p2(2)], ...
        'LineWidth', 2.2, 'LineStyle','-', ...
        'Color', cols(mod(k-1,size(cols,1))+1,:));
end

hLegNode = plot(nan,nan,'o','MarkerFaceColor','k','MarkerEdgeColor','w','LineStyle','none');
hLegSrc  = plot(nan,nan,'p','MarkerFaceColor','r','MarkerEdgeColor','r','LineStyle','none');
hLegLOS  = plot(nan,nan,'-','LineWidth',2.2,'Color','k');

legend([hLegNode,hLegSrc,hLegLOS], ...
    {'巡检节点','点源','LOS直线边'}, ...
    'Location','northeastoutside');
hold off;

% 6.4 剂量分解条形图
[edgeDoseList, edgeTypeList] = getDoseBreakdown(gBestPath, doseMat, edgeType); %#ok<ASGLU>
figure('Name','剂量分解条形图（无障碍-精细模型）','Position',[120,120,1100,420]);
bar(edgeDoseList);
grid on;
xlabel('回路边序号');
ylabel('该边剂量贡献');
title('最优巡检回路的剂量分解（逐边贡献，全部为LOS边）');
legend({'边剂量'}, 'Location','northeast');

%% ===================== 7) Save MAT result =====================
resultFile = 'result_case3_LOS_refined.mat';

save(resultFile, ...
    'gBestPath', 'gBestFitness', 'convergence', ...
    'doseMat', 'lenMat', 'costMat', 'edgeType', ...
    'xCoords', 'yCoords', 'nodeData', 'numNodes', ...
    'sources', 'sourceStrength', ...
    'epsLen', 'dsLine', ...
    'xMin', 'xMax', 'yMin', 'yMax', ...
    'xVec', 'yVec', 'XX', 'YY', ...
    'radGrid', ...
    'M', 'maxIter', 'eliteKeep', ...
    'twoOptIters_defense', 'twoOptIters_attack', 'twoOptIters_escape', ...
    'popBall', 'popBreed', 'popBaby', ...
    'defenseProb0', 'attackProb0', 'escapeProb0', ...
    'totalDose', 'totalLen', 'nLOSEdges', ...
    'edgeDoseList', 'edgeTypeList');

fprintf('\n已保存 MAT 结果文件: %s\n', resultFile);
%% ============================================================
%                      Functions
%% ============================================================

function c = tourCost(path, mat)
    c = 0;
    for i = 1:(length(path)-1)
        c = c + mat(path(i), path(i+1));
    end
end

function [doseList, typeList] = getDoseBreakdown(tour, doseMat, edgeType)
    m = numel(tour)-1;
    doseList = zeros(m,1);
    typeList = zeros(m,1);
    for k = 1:m
        a = tour(k); b = tour(k+1);
        doseList(k) = doseMat(a,b);
        typeList(k) = edgeType(a,b);
    end
end

%% ---------- Straight-line dose integral ----------
function dose = lineDoseIntegral(p1, p2, sources, strength, ds)
    if nargin < 5, ds = 1.0; end
    L = hypot(p2(1)-p1(1), p2(2)-p1(2));
    n = max(2, ceil(L/ds)+1);
    t = linspace(0,1,n);
    x = p1(1) + t*(p2(1)-p1(1));
    y = p1(2) + t*(p2(2)-p1(2));

    rad = zeros(size(x));
    for s = 1:size(sources,1)
        dx = x - sources(s,1);
        dy = y - sources(s,2);
        rad = rad + strength(s) ./ (dx.^2 + dy.^2 + 10);
    end

    ds_eff = L/(n-1);
    dose = trapz(rad) * ds_eff / 100;
end

%% ---------- CPO behaviors & permutation ops ----------
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

    if rand() < pSwap
        p = swapPerm(p);
    end
    if rand() < pInsert
        p = insertPerm(p);
    end
    if rand() < pInv
        p = inversionPerm(p);
    end

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

function path2 = twoOptLocal(path, costMat, iters)
    if iters <= 0
        path2 = path;
        return;
    end
    p = path;
    n = numel(p);
    for k = 1:iters
        i = randi([2, n-3]);
        j = randi([i+1, n-2]);

        old = costMat(p(i-1), p(i)) + costMat(p(j), p(j+1));
        neu = costMat(p(i-1), p(j)) + costMat(p(i), p(j+1));

        if neu < old
            p(i:j) = fliplr(p(i:j));
        end
    end
    path2 = p;
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
