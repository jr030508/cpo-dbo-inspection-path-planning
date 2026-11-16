%% ============================================================
%  Case3 (51 nodes) Nuclear scenario dose-aware routing
%  DBO (Dung Beetle Optimizer) - PERMUTATION-SAFE version
%
%  What is fixed vs your original Case3 DBO:
%   1) PERMUTATION-SAFE operators for TSP-like tour
%      - No segment overwrite that creates duplicates/missing nodes
%      - Use swap/insert/inversion/scramble mutations (all keep permutation)
%   2) Deterministic repair (NO random reset)
%   3) Keep your DBO 3-subpopulation structure:
%      - Ball beetles: exploitation (mild mutation + 2-opt)
%      - Breeding beetles: guided exploration (OX crossover to gBest + mutation + 2-opt)
%      - Baby beetles: diversification (strong mutation + optional 2-opt)
%   4) Light 2-opt local improvement for stability
%   5) Radiation heatmap uses SAME sourceStrength & decay as dose model
%
%  Tour definition: start at node 1, visit all nodes once, return to 1.
%% ============================================================

clc;
clear;
close all;

rng(7); % for stability test use rng('shuffle')

%% ===================== 51 nodes (Case3) =====================
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

nodeIDs = nodeData(:,1);
xCoords = nodeData(:,2);
yCoords = nodeData(:,3);
numNodes = length(nodeIDs);

%% ===================== 7 sources (Case3) ====================
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

%% ===================== Build doseMatrix (μSv) =================
doseMatrix = zeros(numNodes, numNodes);

nodeRadLevel = zeros(numNodes, 1);
for i = 1:numNodes
    totalRad = 0;
    for s = 1:size(sources, 1)
        dist = norm([xCoords(i)-sources(s,1), yCoords(i)-sources(s,2)]);
        totalRad = totalRad + sourceStrength(s) / (dist^2 + 10);
    end
    nodeRadLevel(i) = totalRad;
end

for i = 1:numNodes
    for j = i+1:numNodes
        dist = norm([xCoords(i)-xCoords(j), yCoords(i)-yCoords(j)]);
        avgRad = (nodeRadLevel(i) + nodeRadLevel(j)) / 2;

        midX = (xCoords(i) + xCoords(j)) / 2;
        midY = (yCoords(i) + yCoords(j)) / 2;
        midRad = 0;
        for s = 1:size(sources, 1)
            distToSource = norm([midX-sources(s,1), midY-sources(s,2)]);
            midRad = midRad + (sourceStrength(s)/2) / (distToSource^2 + 10);
        end

        dose = dist * (avgRad + midRad) / 100;
        doseMatrix(i,j) = dose;
        doseMatrix(j,i) = dose;
    end
end

% hotspot edges: top 5% -> x3
allDoses = doseMatrix(:);
allDoses(allDoses == 0) = [];
threshold = prctile(allDoses, 95);

hotSpots = [];
for i = 1:numNodes
    for j = i+1:numNodes
        if doseMatrix(i,j) > threshold
            hotSpots = [hotSpots; i, j]; %#ok<AGROW>
            doseMatrix(i,j) = doseMatrix(i,j) * 3;
            doseMatrix(j,i) = doseMatrix(j,i) * 3;
        end
    end
end

doseMatrix(logical(eye(numNodes))) = 0;

%% ===== Print doseMatrix (optional, big) =====
fprintf('\n===== 节点间有效剂量矩阵 (单位: μSv) =====\n');
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

%% ===================== DBO parameters ========================
M = 400;
maxIter = 5000;
eliteKeep = max(2, round(0.02*M)); % keep 2% elites for stability

popBall  = round(M * 0.50);
popBreed = round(M * 0.25);
popBaby  = M - popBall - popBreed;

% local improvement intensity (small is enough)
twoOpt_ball  = 0;
twoOpt_breed = 0;
twoOpt_baby  = 0;

%% ===================== Initialize population =================
population = cell(M, 1);
fitness    = zeros(M, 1);

for i = 1:M
    perm = randperm(numNodes-1) + 1;     % 2..numNodes
    population{i} = [1, perm, 1];
    fitness(i) = calculateFitness(population{i}, doseMatrix);
end

pBestPath = population;
pBestFitness = fitness;
[gBestFitness, gBestIdx] = min(fitness);
gBestPath = population{gBestIdx};

convergence = zeros(maxIter, 1);
convergence(1) = gBestFitness;

%% ===================== Main loop (DBO) =======================
fprintf('\n开始蜣螂优化(DBO)算法（Permutation-safe, Case3 51节点）...\n');

for iter = 1:maxIter
    t = iter / maxIter;

    % (Optional) dynamic factors to gradually shift exploration->exploitation
    exploreScale = max(0.10, 1.00 - 0.90*t);  % decreases over time
    exploitScale = min(1.00, 0.30 + 0.70*t);  % increases over time

    % elites (global stability)
    [~, idxSort] = sort(fitness, 'ascend');
    eliteIdx = idxSort(1:eliteKeep);

    newPopulation = population;
    newFitness = fitness;
    for e = 1:eliteKeep
        newPopulation{e} = population{eliteIdx(e)};
        newFitness(e) = fitness(eliteIdx(e));
    end

    % 1) Ball beetles: exploitation (mild mutation around pBest)
    for i = eliteKeep+1:popBall
        current = population{i};
        currentFit = fitness(i);

        base = pBestPath{i};

        % mild, exploitation-oriented mutations
        cand = base;
        cand = mutate_perm(cand, numNodes, ...
            0.40*exploitScale, 0.10*exploitScale, 0.15*exploitScale, 0.05*exploitScale);
        cand = twoOptLocal(cand, doseMatrix, twoOpt_ball);

        cand = validatePath_det(cand, numNodes);
        candFit = calculateFitness(cand, doseMatrix);

        if candFit < currentFit
            newPopulation{i} = cand;
            newFitness(i) = candFit;

            if candFit < pBestFitness(i)
                pBestPath{i} = cand;
                pBestFitness(i) = candFit;
            end
        else
            newPopulation{i} = current;
            newFitness(i) = currentFit;
        end
    end

    % 2) Breeding beetles: guided exploration (OX crossover to gBest + mutation)
    for i = max(popBall+1, eliteKeep+1):(popBall+popBreed)
        current = population{i};
        currentFit = fitness(i);

        base = pBestPath{i};

        % guided crossover with gBest (permutation-safe)
        if rand() < 0.65*exploreScale + 0.15
            cand = oxCrossover(base, gBestPath);
        else
            cand = base;
        end

        % moderate mutations
        cand = mutate_perm(cand, numNodes, ...
            0.25*exploreScale, 0.20*exploreScale, 0.25*exploreScale, 0.10*exploreScale);
        cand = twoOptLocal(cand, doseMatrix, twoOpt_breed);

        cand = validatePath_det(cand, numNodes);
        candFit = calculateFitness(cand, doseMatrix);

        if candFit < currentFit
            newPopulation{i} = cand;
            newFitness(i) = candFit;

            if candFit < pBestFitness(i)
                pBestPath{i} = cand;
                pBestFitness(i) = candFit;
            end
        else
            newPopulation{i} = current;
            newFitness(i) = currentFit;
        end
    end

    % 3) Baby beetles: diversification (strong perturbations)
    for i = max(popBall+popBreed+1, eliteKeep+1):M
        current = population{i};
        currentFit = fitness(i);

        base = pBestPath{i};

        cand = base;

        % strong diversify: scramble + inversion more likely
        cand = mutate_perm(cand, numNodes, ...
            0.20*exploreScale, 0.25*exploreScale, 0.45*exploreScale, 0.35*exploreScale);

        % sometimes restart partially from gBest but via OX (still safe)
        if rand() < 0.25*exploreScale
            cand = oxCrossover(cand, gBestPath);
        end

        cand = twoOptLocal(cand, doseMatrix, twoOpt_baby);

        cand = validatePath_det(cand, numNodes);
        candFit = calculateFitness(cand, doseMatrix);

        if candFit < currentFit
            newPopulation{i} = cand;
            newFitness(i) = candFit;

            if candFit < pBestFitness(i)
                pBestPath{i} = cand;
                pBestFitness(i) = candFit;
            end
        else
            newPopulation{i} = current;
            newFitness(i) = currentFit;
        end
    end

    population = newPopulation;
    fitness = newFitness;

    % update gBest
    [minFit, minIdx] = min(fitness);
    if minFit < gBestFitness
        gBestFitness = minFit;
        gBestPath = population{minIdx};
    end

    convergence(iter) = gBestFitness;

    if mod(iter, 50) == 0
        fprintf('DBO迭代 %d/%d | Best dose: %.4f μSv\n', iter, maxIter, gBestFitness);
    end
end

%% ===================== Visualization =========================

% convergence curve
figure('Name', '算法收敛曲线', 'Position', [100, 100, 800, 400]);
plot(convergence, 'LineWidth', 2);
xlabel('迭代次数');
ylabel('最优剂量 (μSv)');
title('蜣螂优化器(DBO)收敛曲线 (Case3: 51节点, Permutation-safe)');
grid on;
legend('最优剂量', 'Location', 'northeast');

% nodes & sources
figure('Name', '节点与放射源分布', 'Position', [100, 100, 900, 700]);
hold on; grid on; axis equal;
xlim([0 200]); ylim([0 200]);
set(gca, 'FontSize', 12);
xlabel('X坐标 (m)', 'FontSize', 14);
ylabel('Y坐标 (m)', 'FontSize', 14);
title('节点与放射源分布图 (案例3: 51节点)', 'FontSize', 16);

scatter(xCoords, yCoords, 100, 'filled', 'MarkerFaceColor', [0.1, 0.5, 0.8], ...
    'MarkerEdgeColor', [0, 0.2, 0.4], 'LineWidth', 1.5);

offset = 3;
for i = 1:numNodes
    textPos = [xCoords(i) + offset, yCoords(i) + offset/2];
    text(textPos(1), textPos(2), sprintf('%d', nodeIDs(i)), ...
        'FontSize', 10, 'FontWeight', 'bold', 'Color', 'k', ...
        'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle');
end

scatter(sources(:,1), sources(:,2), 200, 'r', 'pentagram', 'filled', 'LineWidth', 1.5);
for s = 1:size(sources, 1)
    text(sources(s,1), sources(s,2)-5, sprintf('放射源%d', s),...
        'FontSize', 12, 'FontWeight', 'bold', 'Color', 'r', ...
        'HorizontalAlignment', 'center', 'VerticalAlignment', 'top');
end
legend('节点位置', '放射源位置', 'Location', 'northeast', 'FontSize', 12);
box on; grid minor;
set(gca, 'GridColor', [0.85, 0.85, 0.85], 'GridAlpha', 0.7);
hold off;

% best tour plot
figure('Name', '最优路径规划', 'Position', [100, 100, 900, 700]);
hold on; grid on; axis equal;
xlim([0 200]); ylim([0 200]);
set(gca, 'FontSize', 12);
xlabel('X坐标 (m)', 'FontSize', 14);
ylabel('Y坐标 (m)', 'FontSize', 14);
title(sprintf('案例3最优路径规划 (总剂量: %.4f μSv)', gBestFitness), 'FontSize', 16);

scatter(xCoords, yCoords, 100, 'filled', 'MarkerFaceColor', [0.1, 0.5, 0.8], ...
    'MarkerEdgeColor', [0, 0.2, 0.4], 'LineWidth', 1.5);

for i = 1:numNodes
    textPos = [xCoords(i) + offset, yCoords(i) + offset/2];
    text(textPos(1), textPos(2), sprintf('%d', nodeIDs(i)), ...
        'FontSize', 10, 'FontWeight', 'bold', 'Color', 'k', ...
        'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle');
end

scatter(sources(:,1), sources(:,2), 200, 'r', 'pentagram', 'filled', 'LineWidth', 1.5);

pathColors = jet(length(gBestPath)-1);
for k = 1:(length(gBestPath)-1)
    a = gBestPath(k);
    b = gBestPath(k+1);
    plot([xCoords(a), xCoords(b)], [yCoords(a), yCoords(b)], ...
         'Color', pathColors(k,:), 'LineWidth', 2);
    if mod(k, 5) == 0
        midX = (xCoords(a) + xCoords(b)) / 2;
        midY = (yCoords(a) + yCoords(b)) / 2;
        text(midX, midY, sprintf('%d→%d', a, b), ...
            'FontSize', 8, 'BackgroundColor', [1,1,1,0.7], ...
            'HorizontalAlignment', 'center');
    end
end

plot(xCoords(1), yCoords(1), 'o', 'MarkerSize', 12, 'MarkerFaceColor', 'g', 'MarkerEdgeColor', 'k');
text(xCoords(1), yCoords(1)+5, '起点', 'FontSize', 10, 'FontWeight', 'bold', 'Color', 'g', ...
    'HorizontalAlignment', 'center');

legend('节点位置', '放射源位置', '路径', 'Location', 'northeast', 'FontSize', 12);
box on; grid minor;
set(gca, 'GridColor', [0.85, 0.85, 0.85], 'GridAlpha', 0.7);
hold off;

% dose heatmap
figure('Name', '节点间有效剂量热力图', 'Position', [100, 100, 1000, 850]);

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
c.Label.String = '有效剂量 (μSv)';
c.Label.FontSize = 12;
c.FontSize = 11;

set(gca, 'XTick', 1:numNodes, 'YTick', 1:numNodes, ...
         'FontSize', 8, 'XAxisLocation', 'top');
xlabel('目标节点', 'FontSize', 14, 'FontWeight', 'bold');
ylabel('起始节点', 'FontSize', 14, 'FontWeight', 'bold');
title('节点间有效剂量热力图 (案例3: 51节点)', 'FontSize', 16, 'FontWeight', 'bold');

hold on;
for i = 1:numNodes+1
    plot([0.5, numNodes+0.5], [i-0.5, i-0.5], 'k-', 'LineWidth', 0.5);
    plot([i-0.5, i-0.5], [0.5, numNodes+0.5], 'k-', 'LineWidth', 0.5);
end

p95 = prctile(allDoses, 95);
doseValues = doseMatrix;
doseValues(doseValues == 0) = NaN;

for i = 1:numNodes
    for j = 1:numNodes
        if ~isnan(doseValues(i,j)) && doseValues(i,j) > p95
            text(j, i, sprintf('%.0f', doseValues(i,j)), ...
                 'HorizontalAlignment', 'center', ...
                 'FontSize', 6, 'FontWeight', 'bold', 'Color', 'k');
        end
    end
end

for k = 1:size(hotSpots, 1)
    nodeA = hotSpots(k, 1);
    nodeB = hotSpots(k, 2);
    rectangle('Position', [nodeB-0.5, nodeA-0.5, 1, 1], ...
              'EdgeColor', 'r', 'LineWidth', 1.5, 'LineStyle', '--');
end

annotation('textbox', [0.15, 0.02, 0.7, 0.05], ...
           'String', '注: 红色虚线框表示高辐射路径（前5%），数值仅显示 >95% 分位数的剂量', ...
           'FontSize', 11, 'EdgeColor', 'none', 'HorizontalAlignment', 'center');

axis square;
set(gcf, 'Color', 'w');
hold off;

% radiation field heatmap (CONSISTENT)
figure('Name', '辐射环境热力图', 'Position', [100, 100, 900, 700]);
hold on;

[xGrid, yGrid] = meshgrid(0:2:200, 0:2:200);
radGrid = zeros(size(xGrid));

for i = 1:size(xGrid,1)
    for j = 1:size(xGrid,2)
        totalRad = 0;
        for s = 1:size(sources, 1)
            dist = norm([xGrid(i,j)-sources(s,1), yGrid(i,j)-sources(s,2)]);
            totalRad = totalRad + sourceStrength(s) / (dist^2 + 10);
        end
        radGrid(i,j) = totalRad;
    end
end

contourf(xGrid, yGrid, radGrid, 30, 'LineStyle', 'none');
colormap(jet);
c = colorbar;
c.Label.String = '辐射水平';
c.Label.FontSize = 12;
title('辐射环境热力图 (案例3: 51节点)', 'FontSize', 16);
xlabel('X坐标 (m)', 'FontSize', 14);
ylabel('Y坐标 (m)', 'FontSize', 14);

scatter(xCoords, yCoords, 80, 'k', 'filled', 'MarkerEdgeColor', 'w', 'LineWidth', 1.5);
for i = 1:numNodes
    text(xCoords(i)+3, yCoords(i)+3, sprintf('%d', nodeIDs(i)), ...
         'FontSize', 9, 'FontWeight', 'bold', 'Color', 'w', ...
         'BackgroundColor', [0,0,0,0.5], ...
         'HorizontalAlignment', 'center');
end

scatter(sources(:,1), sources(:,2), 150, 'r', 'pentagram', 'filled', 'LineWidth', 1.5);
contour(xGrid, yGrid, radGrid, 10, 'LineColor', 'k', 'LineWidth', 0.5);

legend('节点位置', '放射源位置', 'Location', 'northeast', 'FontSize', 12);
axis equal;
xlim([0 200]); ylim([0 200]);
hold off;

%% ===================== Output summary ========================
fprintf('\n===== 案例3复现结果 (DBO Permutation-safe) =====\n');
fprintf('最优路径: %s\n', strjoin(cellstr(num2str(gBestPath')), '->'));
fprintf('最优路径总剂量: %.4f μSv\n', gBestFitness);
fprintf('算法总迭代次数: %d\n', maxIter);

fprintf('放射源位置:\n');
for s = 1:size(sources, 1)
    fprintf('放射源%d: (%.1f, %.1f), Strength=%.0f\n', s, sources(s,1), sources(s,2), sourceStrength(s));
end

%% ===================== Functions ============================

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

    i = randi([1, n-1]);
    j = randi([i+1, n]);

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

% fitness: total dose
function totalDose = calculateFitness(path, doseMatrix)
    totalDose = 0;
    for i = 1:(length(path)-1)
        totalDose = totalDose + doseMatrix(path(i), path(i+1));
    end
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