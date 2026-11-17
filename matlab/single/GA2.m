%% ============================================================
%  Case3 (51 nodes) Nuclear scenario dose-aware routing
%  GA (Genetic Algorithm) - PERMUTATION-SAFE version
%
%  Fixes vs your original Case3 GA:
%   1) Initialization is always a valid tour (no redundant "补缺" that breaks length)
%   2) OX crossover kept (permutation-safe)
%   3) Mutation uses permutation-safe operators (swap / inversion)
%   4) Deterministic repair validatePath_det() (NO random reset)
%   5) Optional light 2-opt applied to elites (stabilizes convergence)
%   6) Radiation heatmap uses SAME sourceStrength & decay as dose model
%
%  Tour definition: start at node 1, visit all nodes once, return to 1.
%% ============================================================

clc;
clear;
close all;

rng(7);  % use rng('shuffle') to test robustness

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

%% ===================== Build doseMatrix (a.u.) =================
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

%% ===================== GA parameters =========================
popSize = 400;
maxGen  = 5000;
pc      = 0.85;
pm      = 0.15;

eliteRatio = 0.02;
eliteNum   = max(2, round(popSize * eliteRatio));

% optional: apply 2-opt to elites (0/1 small is enough)
twoOptEliteIters = 1;

%% ===================== Initialize population =================
population = cell(popSize, 1);
fitness = zeros(popSize, 1);

for i = 1:popSize
    population{i} = initializePath_det(numNodes);
    fitness(i) = calculateFitness(population{i}, doseMatrix);
end

[bestFitness, bestIdx] = min(fitness);
bestPath = population{bestIdx};

convergence = zeros(maxGen, 1);

%% ===================== GA main loop ==========================
fprintf('\n开始遗传算法优化（Permutation-safe, Case3 51节点）...\n');

for gen = 1:maxGen
    % sort current population by fitness (ascending)
    [fitnessSorted, idxSort] = sort(fitness, 'ascend');
    populationSorted = population(idxSort);

    % elites
    elites = populationSorted(1:eliteNum);
    eliteFit = fitnessSorted(1:eliteNum);

    % optionally improve elites with 2-opt (more stable)
    if twoOptEliteIters > 0
        for e = 1:eliteNum
            elites{e} = twoOptLocal(elites{e}, doseMatrix, twoOptEliteIters);
            eliteFit(e) = calculateFitness(elites{e}, doseMatrix);
        end
        % re-sort elites after 2-opt
        [eliteFit, eIdx] = sort(eliteFit, 'ascend');
        elites = elites(eIdx);
    end

    % selection (tournament) from sorted population
    matingPool = selectionTournament(populationSorted, fitnessSorted, popSize);

    % crossover (OX) -> offspring
    offspring = crossoverOX(matingPool, pc, numNodes);

    % mutation -> offspring
    offspring = mutationPerm(offspring, pm, numNodes);

    % evaluate offspring
    offFit = zeros(popSize, 1);
    for i = 1:popSize
        offFit(i) = calculateFitness(offspring{i}, doseMatrix);
    end

    % combine offspring + elites -> next generation
    allPop = [offspring; elites];
    allFit = [offFit; eliteFit];

    [allFitSorted, allIdx] = sort(allFit, 'ascend');
    population = allPop(allIdx(1:popSize));
    fitness    = allFitSorted(1:popSize);

    % update global best
    if fitness(1) < bestFitness
        bestFitness = fitness(1);
        bestPath = population{1};
    end

    convergence(gen) = bestFitness;

    if mod(gen, 100) == 0
        fprintf('代数 %d/%d, 当前最优暴露代价: %.4f a.u.\n', gen, maxGen, bestFitness);
    end
end

%% ===================== Visualization =========================
% convergence curve
figure('Name', '算法收敛曲线', 'Position', [100, 100, 800, 400]);
plot(convergence, 'LineWidth', 2);
xlabel('迭代次数');
ylabel('最优暴露代价 (a.u.)');
title('遗传算法(GA)收敛曲线 (Case3: 51节点, Permutation-safe)');
grid on;
legend('最优暴露代价', 'Location', 'northeast');

% nodes & sources
figure('Name', '节点与放射源分布', 'Position', [100, 100, 900, 700]);
hold on; grid on; axis equal;
xlim([0 200]); ylim([0 200]);
set(gca, 'FontSize', 12);
xlabel('X坐标 (m)', 'FontSize', 14);
ylabel('Y坐标 (m)', 'FontSize', 14);
title('节点与放射源分布图 (案例3)', 'FontSize', 16);

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
for i = 1:size(sources, 1)
    text(sources(i,1), sources(i,2)-5, sprintf('放射源%d', i),...
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
title(sprintf('案例3最优路径规划 (总暴露代价: %.4f a.u.)', bestFitness), 'FontSize', 16);

scatter(xCoords, yCoords, 100, 'filled', 'MarkerFaceColor', [0.1, 0.5, 0.8], ...
    'MarkerEdgeColor', [0, 0.2, 0.4], 'LineWidth', 1.5);

for i = 1:numNodes
    textPos = [xCoords(i) + offset, yCoords(i) + offset/2];
    text(textPos(1), textPos(2), sprintf('%d', nodeIDs(i)), ...
        'FontSize', 10, 'FontWeight', 'bold', 'Color', 'k', ...
        'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle');
end

scatter(sources(:,1), sources(:,2), 200, 'r', 'pentagram', 'filled', 'LineWidth', 1.5);

pathColors = jet(length(bestPath)-1);
for k = 1:(length(bestPath)-1)
    a = bestPath(k);
    b = bestPath(k+1);
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
text(xCoords(1), yCoords(1)+5, '起点', 'FontSize', 12, 'FontWeight', 'bold', 'Color', 'g', ...
    'HorizontalAlignment', 'center');

legend('节点位置', '放射源位置', '路径', 'Location', 'northeast', 'FontSize', 12);
box on; grid minor;
set(gca, 'GridColor', [0.85, 0.85, 0.85], 'GridAlpha', 0.7);
hold off;

% dose heatmap
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
         'FontSize', 8, 'XAxisLocation', 'top');
xlabel('目标节点', 'FontSize', 14, 'FontWeight', 'bold');
ylabel('起始节点', 'FontSize', 14, 'FontWeight', 'bold');
title('节点间暴露代价矩阵热力图 (案例3)', 'FontSize', 16, 'FontWeight', 'bold');

hold on;
for i = 1:numNodes+1
    plot([0.5, numNodes+0.5], [i-0.5, i-0.5], 'k-', 'LineWidth', 0.5);
    plot([i-0.5, i-0.5], [0.5, numNodes+0.5], 'k-', 'LineWidth', 0.5);
end

for i = 1:size(hotSpots, 1)
    nodeA = hotSpots(i, 1);
    nodeB = hotSpots(i, 2);
    rectangle('Position', [nodeB-0.5, nodeA-0.5, 1, 1], ...
              'EdgeColor', 'r', 'LineWidth', 1, 'LineStyle', '--');
end

annotation('textbox', [0.15, 0.02, 0.7, 0.05], ...
           'String', '注: 红色虚线框表示高辐射路径（前5%）', ...
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
fprintf('\n===== 案例3遗传算法结果 (GA Permutation-safe) =====\n');
fprintf('最优路径: %s\n', strjoin(cellstr(num2str(bestPath')), '->'));
fprintf('最优路径总暴露代价: %.4f a.u.\n', bestFitness);
fprintf('算法总迭代次数: %d\n', maxGen);
fprintf('放射源位置:\n');
for i = 1:size(sources, 1)
    fprintf('放射源%d: (%.1f, %.1f), Strength=%.0f\n', i, sources(i,1), sources(i,2), sourceStrength(i));
end

%% ===================== GA helper functions ===================

% Deterministic initialization (always valid)
function path = initializePath_det(numNodes)
    path = [1, randperm(numNodes-1)+1, 1];
end

% Fitness: total dose
function totalDose = calculateFitness(path, doseMatrix)
    totalDose = 0;
    for i = 1:(length(path)-1)
        totalDose = totalDose + doseMatrix(path(i), path(i+1));
    end
end

% Tournament selection (works on sorted or unsorted)
function newPop = selectionTournament(population, fitness, popSize)
    newPop = cell(popSize, 1);
    tournamentSize = 3;
    n = numel(population);

    for i = 1:popSize
        cand = randperm(n, tournamentSize);
        [~, bestLocal] = min(fitness(cand));
        newPop{i} = population{cand(bestLocal)};
    end
end

% OX crossover over full path (keeps 1 at ends)
function offspring = crossoverOX(population, pc, numNodes)
    offspring = population;
    n = numel(population);

    % pairwise
    for i = 1:2:n-1
        if rand < pc
            p1 = population{i};
            p2 = population{i+1};

            a = p1(2:end-1);
            b = p2(2:end-1);

            c1 = oxPerm(a, b);
            c2 = oxPerm(b, a);

            child1 = [1, c1, 1];
            child2 = [1, c2, 1];

            offspring{i}   = validatePath_det(child1, numNodes);
            offspring{i+1} = validatePath_det(child2, numNodes);
        end
    end
end

% OX on permutation only
function child = oxPerm(a, b)
    n = numel(a);
    child = zeros(1, n);

    pts = sort(randperm(n, 2));
    i = pts(1); j = pts(2);

    child(i:j) = a(i:j);

    fillPos = [1:i-1, j+1:n];
    bSeq = b(~ismember(b, child(i:j)));

    child(fillPos) = bSeq;
end

% Mutation: swap or inversion (both permutation-safe)
function pop2 = mutationPerm(population, pm, numNodes)
    pop2 = population;
    n = numel(population);

    for i = 1:n
        if rand < pm
            path = population{i};
            mid = path(2:end-1);

            if rand < 0.50
                % swap
                idx = randperm(numel(mid), 2);
                tmp = mid(idx(1));
                mid(idx(1)) = mid(idx(2));
                mid(idx(2)) = tmp;
            else
                % inversion
                pts = sort(randperm(numel(mid), 2));
                mid(pts(1):pts(2)) = fliplr(mid(pts(1):pts(2)));
            end

            newPath = [1, mid, 1];
            pop2{i} = validatePath_det(newPath, numNodes);
        end
    end
end

% Deterministic repair (no random reset)
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

% Light 2-opt local improvement
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