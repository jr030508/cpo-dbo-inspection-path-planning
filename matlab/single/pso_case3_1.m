%% ============================================================
%  Case2 (30 nodes) Nuclear scenario dose-aware routing
%  IPSO (Chaotic initialization + Discrete PSO-like operators)
%  PERMUTATION-SAFE version (ported from your 51-node final IPSO)
%
%  Keep unchanged:
%   - 30-node coordinates, 4 sources, dose model, manual hotspots, all plots
%  Replace:
%   - CPO-DBO main loop -> IPSO chaos phase + PSO-like phase
%
%  Tour: start at node 1, visit all nodes once, return to 1.
%% ============================================================

%% 清理环境
clc;
clear;
close all;

rng(7);   % 需要看多次运行稳定性可改成 rng('shuffle')

%% 节点坐标数据（表4） 30 nodes
nodeData = [
    1, 60, 0;
    2, 68, 15;
    3, 73, 29;
    4, 52, 44;
    5, 64, 64;
    6, 44, 68;
    7, 24, 76;
    8, 8, 63;
    9, 22, 52;
    10, 32, 48;
    11, 12, 36;
    12, 23, 29;
    13, 12, 12;
    14, 25, 4;
    15, 41, 32;
    16, 30, 65;
    17, 44, 22;
    18, 37, 15;
    19, 43, 57;
    20, 72, 77;
    21, 78, 49;
    22, 26, 41;
    23, 49, 18;
    24, 16, 46;
    25, 48, 60;
    26, 66, 40;
    27, 53, 65;
    28, 57, 32;
    29, 23, 40;
    30, 15, 20;
];

nodeIDs = nodeData(:,1);
xCoords = nodeData(:,2);
yCoords = nodeData(:,3);
numNodes = length(nodeIDs);

%% 假设放射源位置（基于图5） 4 sources
sources = [
    20, 70;  % 放射源1 (左上)
    40, 40;  % 放射源2 (中心)
    60, 60;  % 放射源3 (右上)
    70, 20;  % 放射源4 (右下)
];

%% ===================== 计算剂量矩阵 ==========================
doseMatrix = zeros(numNodes, numNodes);

% 每个节点的辐射水平
nodeRadLevel = zeros(numNodes, 1);
for i = 1:numNodes
    totalRad = 0;
    for s = 1:size(sources, 1)
        dist = norm([xCoords(i)-sources(s,1), yCoords(i)-sources(s,2)]);
        totalRad = totalRad + 1000 / (dist^2 + 5);
    end
    nodeRadLevel(i) = totalRad;
end

% 节点间边剂量
for i = 1:numNodes
    for j = i+1:numNodes
        dist = norm([xCoords(i)-xCoords(j), yCoords(i)-yCoords(j)]);
        avgRad = (nodeRadLevel(i) + nodeRadLevel(j)) / 2;

        midX = (xCoords(i) + xCoords(j)) / 2;
        midY = (yCoords(i) + yCoords(j)) / 2;
        midRad = 0;
        for s = 1:size(sources, 1)
            distToSource = norm([midX-sources(s,1), midY-sources(s,2)]);
            midRad = midRad + 500 / (distToSource^2 + 5);
        end

        dose = dist * (avgRad + midRad) / 100;
        doseMatrix(i,j) = dose;
        doseMatrix(j,i) = dose;
    end
end

% 添加高辐射热点（基于图5）——保持你原来的手工热点
hotSpots = [6, 20; 7, 20; 20, 26; 26, 21];
for k = 1:size(hotSpots, 1)
    a = hotSpots(k, 1);
    b = hotSpots(k, 2);
    doseMatrix(a,b) = doseMatrix(a,b) * 5;
    doseMatrix(b,a) = doseMatrix(b,a) * 5;
end

doseMatrix(logical(eye(numNodes))) = 0;

%% ===== 输出完整剂量矩阵 =====
fprintf('\n===== 节点间暴露代价矩阵 (单位: a.u.) =====\n');
fprintf('       ');
for j = 1:numNodes
    fprintf('%6d', j);
end
fprintf('\n');

for i = 1:numNodes
    fprintf('%4d | ', i);
    for j = 1:numNodes
        if i == j
            fprintf('%6s', '-');
        else
            fprintf('%6.1f', doseMatrix(i,j));
        end
    end
    fprintf('\n');
end

csvwrite('dose_matrix_case2.csv', doseMatrix);
fprintf('\n暴露代价矩阵已保存到文件: dose_matrix_case2.csv\n');

%% ===================== IPSO 参数（替换 CPO-DBO） ==========================
M = 400;

% 仍然保持总迭代次数 = maxIter（不改变你出图形式）
maxIter      = 5000;      % 你原来的 maxIter
maxChaosIter = 200;       % 混沌阶段
maxPSOIter   = maxIter - maxChaosIter;  % PSO-like阶段（保证总长度仍=1000）

% Discrete PSO-like operator probabilities
pF1 = 0.80;   % toward pBest (OX)
pF2 = 0.40;   % toward gBest (OX)
pm  = 0.25;   % mutation probability

% mutation mix probabilities
pSwap  = 0.40;
pInsert= 0.20;
pInv   = 0.25;
pScr   = 0.15;

% optional local improvement
twoOptIters = 1;  % set 0 to disable

%% ===================== Initialize (chaotic keys) =============
% chaos keys correspond to nodes 2..numNodes (length numNodes-1)
chaosKeys = initChaosPopulation(M, numNodes-1);

population = cell(M, 1);
fitness    = zeros(M, 1);

for i = 1:M
    population{i} = chaosKeysToPath(chaosKeys(i, :), numNodes);
    fitness(i) = calculateFitness(population{i}, doseMatrix);
end

pBestPath    = population;
pBestFitness = fitness;

[gBestFitness, gBestIdx] = min(fitness);
gBestPath = population{gBestIdx};

convergence = zeros(maxIter, 1);

%% ===================== Chaos phase ===========================
fprintf('\n开始混沌优化阶段 (Case2-30节点 IPSO)...\n');
for iter = 1:maxChaosIter
    chaosKeys = updateChaos(chaosKeys);

    for i = 1:M
        newPath = chaosKeysToPath(chaosKeys(i, :), numNodes);
        newFit  = calculateFitness(newPath, doseMatrix);

        % update current
        population{i} = newPath;
        fitness(i) = newFit;

        % update pBest
        if newFit < pBestFitness(i)
            pBestPath{i} = newPath;
            pBestFitness(i) = newFit;
        end
    end

    % update gBest
    [minFit, minIdx] = min(pBestFitness);
    if minFit < gBestFitness
        gBestFitness = minFit;
        gBestPath = pBestPath{minIdx};
    end

    convergence(iter) = gBestFitness;

    if mod(iter, 10) == 0
        fprintf('混沌迭代 %d/%d, 当前最优暴露代价: %.4f a.u.\n', iter, maxChaosIter, gBestFitness);
    end
end

%% ===================== PSO-like phase ========================
fprintf('\n开始离散PSO优化阶段 (Case2-30节点 IPSO)...\n');
for iter = 1:maxPSOIter
    t = iter / maxPSOIter;

    % gradually reduce exploration
    pm_now  = max(0.05, pm * (1 - 0.7*t));
    pF1_now = min(0.95, pF1 + 0.10*t);
    pF2_now = min(0.95, pF2 + 0.10*t);

    for i = 1:M
        current = population{i};
        cand = current;

        % F1: toward pBest
        if rand < pF1_now
            cand = oxCrossover(cand, pBestPath{i});
        end

        % F2: toward gBest
        if rand < pF2_now
            cand = oxCrossover(cand, gBestPath);
        end

        % F3: mutation
        if rand < pm_now
            cand = mutate_perm(cand, numNodes, pSwap, pInsert, pInv, pScr);
        end

        % optional local improvement
        if twoOptIters > 0 && rand < 0.30
            cand = twoOptLocal(cand, doseMatrix, twoOptIters);
        end

        cand = validatePath_det(cand, numNodes);
        candFit = calculateFitness(cand, doseMatrix);

        % update current
        population{i} = cand;
        fitness(i) = candFit;

        % update pBest
        if candFit < pBestFitness(i)
            pBestPath{i} = cand;
            pBestFitness(i) = candFit;
        end
    end

    % update gBest
    [minFit, minIdx] = min(pBestFitness);
    if minFit < gBestFitness
        gBestFitness = minFit;
        gBestPath = pBestPath{minIdx};
    end

    convergence(maxChaosIter + iter) = gBestFitness;

    if mod(iter, 50) == 0
        fprintf('PSO迭代 %d/%d, 当前最优暴露代价: %.4f a.u.\n', iter, maxPSOIter, gBestFitness);
    end
end

%% ===================== 结果可视化 ===========================
% 收敛曲线（形式保持不变：单条曲线）
figure('Name', '算法收敛曲线', 'Position', [100, 100, 800, 400]);
plot(convergence, 'LineWidth', 2);
xlabel('迭代次数');
ylabel('最优暴露代价 (a.u.)');
title('混沌初始化 + 离散PSO(IPSO) 收敛曲线 (案例2)');
grid on;
legend('最优暴露代价', 'Location', 'northeast');

%% ===== 增强节点图 =====
figure('Name', '节点与放射源分布', 'Position', [100, 100, 900, 700]);
hold on; grid on; axis equal;
xlim([0 80]); ylim([0 80]);
set(gca, 'FontSize', 12);
xlabel('X坐标 (m)', 'FontSize', 14);
ylabel('Y坐标 (m)', 'FontSize', 14);
title('节点与放射源分布图 (案例2)', 'FontSize', 16);

scatter(xCoords, yCoords, 200, 'filled', 'MarkerFaceColor', [0.1, 0.5, 0.8], ...
    'MarkerEdgeColor', [0, 0.2, 0.4], 'LineWidth', 1.5);

offset = 1.5;
for i = 1:numNodes
    if i == 1 || i == 14
        textPos = [xCoords(i), yCoords(i) - offset*1.5];
    elseif i == 7 || i == 20
        textPos = [xCoords(i) - offset*2, yCoords(i)];
    else
        textPos = [xCoords(i) + offset, yCoords(i) + offset/2];
    end

    text(textPos(1), textPos(2), sprintf('%d', nodeIDs(i)), ...
        'FontSize', 12, 'FontWeight', 'bold', 'Color', 'k', ...
        'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle');
end

scatter(sources(:,1), sources(:,2), 300, 'r', 'pentagram', 'filled', 'LineWidth', 1.5);
text(sources(:,1), sources(:,2)-3, {'放射源1','放射源2','放射源3','放射源4'},...
    'FontSize', 14, 'FontWeight', 'bold', 'Color', 'r', ...
    'HorizontalAlignment', 'center', 'VerticalAlignment', 'top');

legend('节点位置', '放射源位置', 'Location', 'northeast', 'FontSize', 12);
box on; grid minor;
set(gca, 'GridColor', [0.85, 0.85, 0.85], 'GridAlpha', 0.7);
hold off;

%% ===== 最优路径图 =====
figure('Name', '最优路径规划', 'Position', [100, 100, 900, 700]);
hold on; grid on; axis equal;
xlim([0 80]); ylim([0 80]);
set(gca, 'FontSize', 12);
xlabel('X坐标 (m)', 'FontSize', 14);
ylabel('Y坐标 (m)', 'FontSize', 14);
title(sprintf('案例2最优路径规划 (总暴露代价: %.4f a.u.)', gBestFitness), 'FontSize', 16);

scatter(xCoords, yCoords, 200, 'filled', 'MarkerFaceColor', [0.1, 0.5, 0.8], ...
    'MarkerEdgeColor', [0, 0.2, 0.4], 'LineWidth', 1.5);

for i = 1:numNodes
    if i == 1 || i == 14
        textPos = [xCoords(i), yCoords(i) - offset*1.5];
    elseif i == 7 || i == 20
        textPos = [xCoords(i) - offset*2, yCoords(i)];
    else
        textPos = [xCoords(i) + offset, yCoords(i) + offset/2];
    end

    text(textPos(1), textPos(2), sprintf('%d', nodeIDs(i)), ...
        'FontSize', 12, 'FontWeight', 'bold', 'Color', 'k', ...
        'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle');
end

scatter(sources(:,1), sources(:,2), 300, 'r', 'pentagram', 'filled', 'LineWidth', 1.5);
text(sources(:,1), sources(:,2)-3, {'放射源1','放射源2','放射源3','放射源4'},...
    'FontSize', 14, 'FontWeight', 'bold', 'Color', 'r', ...
    'HorizontalAlignment', 'center', 'VerticalAlignment', 'top');

pathColors = jet(length(gBestPath)-1);
for k = 1:(length(gBestPath)-1)
    a = gBestPath(k);
    b = gBestPath(k+1);
    plot([xCoords(a), xCoords(b)], [yCoords(a), yCoords(b)], ...
        'Color', pathColors(k,:), 'LineWidth', 3);
    midX = (xCoords(a) + xCoords(b)) / 2;
    midY = (yCoords(a) + yCoords(b)) / 2;
    text(midX, midY, sprintf('%d→%d', a, b), ...
        'FontSize', 10, 'BackgroundColor', [1,1,1,0.7], ...
        'HorizontalAlignment', 'center');
end

plot(xCoords(1), yCoords(1), 'o', 'MarkerSize', 15, ...
    'MarkerFaceColor', 'g', 'MarkerEdgeColor', 'k');
text(xCoords(1), yCoords(1)+3, '起点', 'FontSize', 12, 'FontWeight', 'bold', ...
    'Color', 'g', 'HorizontalAlignment', 'center');

legend('节点位置', '放射源位置', '路径', 'Location', 'northeast', 'FontSize', 12);
box on; grid minor;
set(gca, 'GridColor', [0.85, 0.85, 0.85], 'GridAlpha', 0.7);
hold off;

%% ===== 增强型剂量热力图 =====
figure('Name', '节点间暴露代价热力图', 'Position', [100, 100, 1000, 850]);

colors = [0.0, 0.0, 1.0;
          0.0, 1.0, 0.0;
          1.0, 1.0, 0.0;
          1.0, 0.5, 0.0;
          1.0, 0.0, 0.0];
colorPos = [0, 0.3, 0.6, 0.8, 1];
customMap = interp1(colorPos, colors, linspace(0,1,256));

imagesc(doseMatrix);
colormap(customMap);
c = colorbar;
c.Label.String = '暴露代价 (a.u.)';
c.Label.FontSize = 12;
c.FontSize = 11;

set(gca, 'XTick', 1:numNodes, 'YTick', 1:numNodes, ...
         'FontSize', 10, 'XAxisLocation', 'top');
xlabel('目标节点', 'FontSize', 14, 'FontWeight', 'bold');
ylabel('起始节点', 'FontSize', 14, 'FontWeight', 'bold');
title('节点间暴露代价矩阵热力图 (案例2)', 'FontSize', 16, 'FontWeight', 'bold');

hold on;
for i = 1:numNodes+1
    plot([0.5, numNodes+0.5], [i-0.5, i-0.5], 'k-', 'LineWidth', 0.5);
    plot([i-0.5, i-0.5], [0.5, numNodes+0.5], 'k-', 'LineWidth', 0.5);
end

doseValues = doseMatrix;
doseValues(doseValues == 0) = NaN;
for i = 1:numNodes
    for j = 1:numNodes
        if ~isnan(doseValues(i,j)) && doseValues(i,j) > 100
            text(j, i, sprintf('%.0f', doseValues(i,j)), ...
                 'HorizontalAlignment', 'center', ...
                 'FontSize', 9, 'FontWeight', 'bold', 'Color', 'k');
        end
    end
end

for k = 1:size(hotSpots, 1)
    nodeA = hotSpots(k, 1);
    nodeB = hotSpots(k, 2);
    rectangle('Position', [nodeB-0.5, nodeA-0.5, 1, 1], ...
              'EdgeColor', 'r', 'LineWidth', 2, 'LineStyle', '--');
end

annotation('textbox', [0.15, 0.02, 0.7, 0.05], ...
           'String', '注: 红色虚线框表示高辐射路径，数值仅显示大于100a.u.的剂量', ...
           'FontSize', 11, 'EdgeColor', 'none', 'HorizontalAlignment', 'center');

axis square;
set(gcf, 'Color', 'w');
hold off;

%% ===== 辐射环境热力图 =====
figure('Name', '辐射环境热力图', 'Position', [100, 100, 900, 700]);
hold on;

[xGrid, yGrid] = meshgrid(0:1:80, 0:1:80);
radGrid = zeros(size(xGrid));

for i = 1:size(xGrid,1)
    for j = 1:size(xGrid,2)
        totalRad = 0;
        for s = 1:size(sources, 1)
            dist = norm([xGrid(i,j)-sources(s,1), yGrid(i,j)-sources(s,2)]);
            totalRad = totalRad + 1000 / (dist^2 + 5);
        end
        radGrid(i,j) = totalRad;
    end
end

contourf(xGrid, yGrid, radGrid, 30, 'LineStyle', 'none');
colormap(jet);
c = colorbar;
c.Label.String = '辐射水平';
c.Label.FontSize = 12;
title('辐射环境热力图 (案例2)', 'FontSize', 16);
xlabel('X坐标 (m)', 'FontSize', 14);
ylabel('Y坐标 (m)', 'FontSize', 14);

scatter(xCoords, yCoords, 100, 'k', 'filled', 'MarkerEdgeColor', 'w', 'LineWidth', 1.5);
for i = 1:numNodes
    text(xCoords(i)+1.5, yCoords(i)+1.5, sprintf('%d', nodeIDs(i)), ...
         'FontSize', 11, 'FontWeight', 'bold', 'Color', 'w', ...
         'BackgroundColor', [0,0,0,0.5], ...
         'HorizontalAlignment', 'center');
end

scatter(sources(:,1), sources(:,2), 200, 'r', 'pentagram', 'filled', 'LineWidth', 1.5);
contour(xGrid, yGrid, radGrid, 10, 'LineColor', 'k', 'LineWidth', 0.5);

legend('节点位置', '放射源位置', 'Location', 'northeast', 'FontSize', 12);
axis equal;
xlim([0 80]); ylim([0 80]);
hold off;

%% 输出结果（格式保持不变）
fprintf('\n===== 案例2复现结果 (Permutation-safe IPSO) =====\n');
fprintf('最优路径: %s\n', strjoin(cellstr(num2str(gBestPath')), '->'));
fprintf('最优路径总暴露代价: %.4f a.u.\n', gBestFitness);
fprintf('算法总迭代次数: %d (混沌:%d + PSO:%d)\n', maxIter, maxChaosIter, maxPSOIter);

%% ===================== 函数区 ==============================

% init chaotic keys in (0,1)
function chaosPop = initChaosPopulation(M, dim)
    chaosPop = zeros(M, dim);
    for i = 1:M
        x = rand(1, dim);
        bad = [0, 0.25, 0.5, 0.75, 1];
        while any(ismembertol(x, bad, 1e-12))
            x = rand(1, dim);
        end
        chaosPop(i, :) = x;
    end
end

% chaos keys -> tour [1, perm(2..N), 1]
function path = chaosKeysToPath(keys, numNodes)
    nodes = 2:numNodes;
    [~, idx] = sort(keys, 'ascend');
    perm = nodes(idx);
    path = [1, perm, 1];
end

% logistic map update
function newChaos = updateChaos(oldChaos)
    newChaos = 4 * oldChaos .* (1 - oldChaos);
end

% fitness: total dose
function totalDose = calculateFitness(path, doseMatrix)
    totalDose = 0;
    for i = 1:(length(path)-1)
        totalDose = totalDose + doseMatrix(path(i), path(i+1));
    end
end

% OX crossover (path-level)
function childPath = oxCrossover(parentA, parentB)
    a = parentA(2:end-1);
    b = parentB(2:end-1);
    c = oxPerm(a, b);
    childPath = [1, c, 1];
end

% OX crossover (perm-level)
function child = oxPerm(a, b)
    n = numel(a);
    child = zeros(1,n);

    pts = sort(randperm(n, 2));
    i = pts(1); j = pts(2);

    child(i:j) = a(i:j);

    fillPos = [1:i-1, j+1:n];
    bSeq = b(~ismember(b, child(i:j)));

    child(fillPos) = bSeq;
end

% mutation pack: swap / insert / inversion / scramble
function path2 = mutate_perm(path, numNodes, pSwap, pInsert, pInv, pScr)
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
    if rand() < pScr
        p = scramblePerm(p);
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
    pts = sort(randperm(n,2));
    p(pts(1):pts(2)) = fliplr(p(pts(1):pts(2)));
end

function p = scramblePerm(p)
    n = numel(p);
    pts = sort(randperm(n,2));
    seg = p(pts(1):pts(2));
    p(pts(1):pts(2)) = seg(randperm(numel(seg)));
end

% deterministic repair (no random reset)
function validPath = validatePath_det(path, numNodes)
    path(1) = 1;
    path(end) = 1;

    mid = path(2:end-1);
    mid = mid(mid>=2 & mid<=numNodes);

    [~, ia] = unique(mid, 'stable');
    mid = mid(sort(ia));

    allNodes = 2:numNodes;
    missing = setdiff(allNodes, mid, 'stable');
    mid = [mid, missing];

    mid = mid(1:(numNodes-1));
    validPath = [1, mid, 1];
end

% light 2-opt local improvement
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

        old = doseMatrix(p(i-1), p(i)) + doseMatrix(p(j), p(j+1));
        new = doseMatrix(p(i-1), p(j)) + doseMatrix(p(i), p(j+1));

        if new < old
            p(i:j) = fliplr(p(i:j));
        end
    end
    path2 = p;
end