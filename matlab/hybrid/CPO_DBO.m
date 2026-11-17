%% ============================================================
%  Case2 (30 nodes) Nuclear scenario dose-aware routing
%  CPO-DBO hybrid (Permutation-safe version)
%
%  Fixes applied (same style as the 51-node fixed version):
%   1) All operators are PERMUTATION-PRESERVING for TSP-like tour
%      - Use OX (Order Crossover) + swap/insert/inversion + scramble
%   2) Dynamic behavior probabilities are ACTUALLY USED
%   3) validate/repair is deterministic and NEVER random-resets a tour
%   4) Add light 2-opt local improvement to stabilize convergence
%   5) Keep your dose model, hotspots, and all visualizations (case2)
%
%  Tour definition: start at node 1, visit all nodes once, return to 1.
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

% 添加高辐射热点（基于图5）
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

%% ===================== CPO-DBO 参数 ==========================
M = 400;
maxIter = 1000;                 % 可先用 5000 做试验
eliteKeep = max(2, round(0.03*M));  % 精英保留 3%

% DBO 分组（保留你原比例，但行为由概率驱动）
popBall  = round(M * 0.50);
popBreed = round(M * 0.25);
popBaby  = M - popBall - popBreed;

% 初始行为概率（随后动态变化并真正参与选择）
defenseProb0 = 0.70;
attackProb0  = 0.25;
escapeProb0  = 0.05;  % 初期少逃跑，后期增加

% 2-opt 局部搜索强度
twoOptIters_defense = 1;
twoOptIters_attack  = 1;
twoOptIters_escape  = 0;

%% ===================== 初始化种群 ===========================
population = cell(M, 1);
fitness    = zeros(M, 1);

for i = 1:M
    perm = randperm(numNodes-1) + 1;     % 2..numNodes 的排列
    population{i} = [1, perm, 1];
    fitness(i) = calculateFitness(population{i}, doseMatrix);
end

pBestPath = population;
pBestFitness = fitness;
[gBestFitness, gBestIdx] = min(fitness);
gBestPath = population{gBestIdx};

convergence = zeros(maxIter, 1);
convergence(1) = gBestFitness;

%% ===================== 主循环 ===============================
fprintf('\n开始 CPO-DBO（Permutation-safe, Case2）...\n');

for iter = 1:maxIter
    t = iter / maxIter;

    % 动态概率（真正用于行为选择）
    defenseProb = defenseProb0 * (1 - 0.6*t);
    attackProb  = attackProb0  * (1 - 0.2*t);
    escapeProb  = min(0.60, escapeProb0 + 0.55*t);

    sProb = defenseProb + attackProb + escapeProb;
    defenseProb = defenseProb / sProb;
    attackProb  = attackProb  / sProb;
    escapeProb  = escapeProb  / sProb;

    % 精英保留
    [~, idxSort] = sort(fitness, 'ascend');
    eliteIdx = idxSort(1:eliteKeep);
    newPopulation = population;
    newFitness = fitness;

    for e = 1:eliteKeep
        newPopulation{e} = population{eliteIdx(e)};
        newFitness(e) = fitness(eliteIdx(e));
    end

    % 更新其余个体
    for i = eliteKeep+1:M
        currentPath = population{i};
        currentFit  = fitness(i);

        % 群体偏置（保持你“滚球/繁殖/小蜣螂”直觉，但不再固定死行为）
        if i <= popBall
            groupBias = 0.10;      % 更偏防御
        elseif i <= popBall + popBreed
            groupBias = -0.05;     % 更偏攻击
        else
            groupBias = -0.10;     % 更偏逃跑
        end

        % 行为选择（动态概率 + 群体偏置）
        r = rand();
        pD = min(0.95, max(0.05, defenseProb + groupBias));
        pA = min(0.95, max(0.05, attackProb  - groupBias/2));
        pS = pD + pA + escapeProb;
        pD = pD / pS; pA = pA / pS; pE = escapeProb / pS;

        if r < pD
            % ---------- 防御：向gBest学习（用OX继承顺序）+ 轻变异 + 2-opt ----------
            newPath = defenseBehavior_perm(currentPath, gBestPath, numNodes);
            newPath = mutate_perm(newPath, numNodes, 0.60, 0.15, 0.15, 0.10);
            newPath = twoOptLocal(newPath, doseMatrix, twoOptIters_defense);
        elseif r < pD + pA
            % ---------- 攻击：OX交叉 + 变异 + 2-opt ----------
            mateIdx = randi(M);
            matePath = pBestPath{mateIdx};
            newPath = oxCrossover(currentPath, matePath);
            newPath = mutate_perm(newPath, numNodes, 0.35, 0.20, 0.30, 0.10);
            newPath = twoOptLocal(newPath, doseMatrix, twoOptIters_attack);
        else
            % ---------- 逃跑：强扰动 + (可选)少量2-opt ----------
            newPath = escapeBehavior_perm(currentPath, numNodes);
            newPath = mutate_perm(newPath, numNodes, 0.25, 0.25, 0.45, 0.25);
            newPath = twoOptLocal(newPath, doseMatrix, twoOptIters_escape);
        end

        % 确保可行（理论上不应破坏，但保底）
        newPath = validatePath_det(newPath, numNodes);
        newFit = calculateFitness(newPath, doseMatrix);

        % 贪婪接受（更稳定）
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

    % 更新全局最优
    [minFitness, minIdx] = min(fitness);
    if minFitness < gBestFitness
        gBestFitness = minFitness;
        gBestPath = population{minIdx};
    end

    convergence(iter) = gBestFitness;

    if mod(iter, 50) == 0
        fprintf('Iter %d/%d | Best dose: %.4f a.u. | (pD=%.2f, pA=%.2f, pE=%.2f)\n',...
            iter, maxIter, gBestFitness, defenseProb, attackProb, escapeProb);
    end
end

%% ===================== 结果可视化 ===========================
% 收敛曲线
figure('Name', '算法收敛曲线', 'Position', [100, 100, 800, 400]);
plot(convergence, 'LineWidth', 2);
xlabel('迭代次数');
ylabel('最优暴露代价 (a.u.)');
title('CPO-DBO（Permutation-safe）收敛曲线 (案例2)');
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

%% 输出结果
fprintf('\n===== 案例2复现结果 (Permutation-safe CPO-DBO) =====\n');
fprintf('最优路径: %s\n', strjoin(cellstr(num2str(gBestPath')), '->'));
fprintf('最优路径总暴露代价: %.4f a.u.\n', gBestFitness);
fprintf('算法总迭代次数: %d\n', maxIter);

%% ===================== 函数区 ==============================

% ---------- 防御（开发）：向gBest学习 + 小扰动 ----------
function newPath = defenseBehavior_perm(path, bestPath, numNodes)
    p = path(2:end-1);
    b = bestPath(2:end-1);

    if rand() < 0.85
        child = oxPerm(p, b);  % 用OX继承顺序信息（可行）
    else
        child = p;
    end

    if rand() < 0.50
        child = swapPerm(child);
    end

    newPath = [1, child, 1];
    newPath = validatePath_det(newPath, numNodes);
end

% ---------- 逃跑（跳出局部）：强扰动 ----------
function newPath = escapeBehavior_perm(path, numNodes)
    p = path(2:end-1);

    if rand() < 0.85
        p = scramblePerm(p);
    end
    if rand() < 0.65
        p = inversionPerm(p);
    end

    newPath = [1, p, 1];
    newPath = validatePath_det(newPath, numNodes);
end

% ---------- OX交叉（Path版本） ----------
function childPath = oxCrossover(parentA, parentB)
    a = parentA(2:end-1);
    b = parentB(2:end-1);
    c = oxPerm(a, b);
    childPath = [1, c, 1];
end

% ---------- OX交叉（Permutation部分） ----------
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

% ---------- 变异组合：swap / insert / inversion / scramble ----------
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

% ---------- 2-opt 局部改进 ----------
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

% ---------- 适应度：路径总剂量 ----------
function totalDose = calculateFitness(path, doseMatrix)
    totalDose = 0;
    for i = 1:(length(path)-1)
        totalDose = totalDose + doseMatrix(path(i), path(i+1));
    end
end

% ---------- 可行性校验（确定性、无随机重置） ----------
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