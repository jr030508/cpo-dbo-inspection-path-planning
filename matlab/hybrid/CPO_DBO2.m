%% ============================================================
%  Case3 (51 nodes) Nuclear scenario dose-aware routing
%  CPO-DBO hybrid (Permutation-safe version)
%
%  Key fixes vs your original code:
%   1) All operators are PERMUTATION-PRESERVING (no duplicate/missing nodes)
%      - Use OX (Order Crossover) + swap/insert/inversion
%   2) "Dynamic probabilities" are ACTUALLY USED to decide behaviors
%      - Each individual chooses defense/attack/escape by current probs
%   3) validate/repair is deterministic and does NOT randomly reset a path
%   4) Add light 2-opt local improvement to stabilize convergence
%
%  NOTE:
%   - This is still a "TSP-like tour": start at node 1, visit all nodes, return 1.
%   - If you want "start->goal shortest path" (not visiting all nodes), tell me,
%     the modeling and encoding are different.
%% ============================================================

%% 清理环境
clc;
clear;
close all;

rng(6);  % 固定随机种子（想看稳定性可改成 rng('shuffle')）

%% 51个节点坐标数据
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

%% 设置7个固定放射源位置（案例3）
sources = [
    35, 35;     % 放射源1
    30, 110;    % 放射源2
    88, 70;     % 放射源3
    88, 142;    % 放射源4
    138, 107;   % 放射源5
    160, 50;    % 放射源6
    168, 180;   % 放射源7
];

% 放射源强度系数
sourceStrength = [8000, 12000, 10000, 15000, 9000, 13000, 11000];

%% 计算剂量矩阵 doseMatrix (μSv)
doseMatrix = zeros(numNodes, numNodes);

% 计算每个节点辐射水平
nodeRadLevel = zeros(numNodes, 1);
for i = 1:numNodes
    totalRad = 0;
    for s = 1:size(sources, 1)
        dist = norm([xCoords(i)-sources(s,1), yCoords(i)-sources(s,2)]);
        totalRad = totalRad + sourceStrength(s) / (dist^2 + 10);
    end
    nodeRadLevel(i) = totalRad;
end

% 计算节点间路径有效剂量
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
doseMatrix(logical(eye(numNodes))) = 0;

% 热点增强（前5%边剂量×3）
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

%% 输出剂量矩阵到CSV
csvwrite('dose_matrix_case3.csv', doseMatrix);
fprintf('\n剂量矩阵已保存到文件: dose_matrix_case3.csv\n');

%% ===================== CPO-DBO参数 ==========================
M = 400;                 % 种群规模
maxIter = 300;         % 迭代次数（可先改小做测试，比如 5000）
eliteKeep = max(2, round(0.02*M));  % 精英保留比例
Tscale = 30000;                 % 固定节奏：以30000为标尺
     

% DBO群体划分（保留你原比例，但行为现在由概率驱动）
popBall  = round(M * 0.4);
popBreed = round(M * 0.2);
popBaby  = M - popBall - popBreed;

% 初始行为概率（随后动态变化并真正参与选择）
defenseProb0 = 0.70;
attackProb0  = 0.25;
escapeProb0  = 0.05;  % 早期少逃跑，后期逐渐增加

% 局部搜索强度（2-opt次数）
twoOptIters_defense = 1;   % 防御更偏开发
twoOptIters_attack  = 1;
twoOptIters_escape  = 0;

%% ===================== 初始化种群 ===========================
population = cell(M,1);
fitness    = zeros(M,1);

for i = 1:M
    perm = randperm(numNodes-1) + 1;         % 2..numNodes 的排列
    population{i} = [1, perm, 1];
    fitness(i) = calculateFitness(population{i}, doseMatrix);
end

pBestPath = population;
pBestFitness = fitness;
[gBestFitness, gBestIdx] = min(fitness);
gBestPath = population{gBestIdx};

convergence = zeros(maxIter, 1);

%% ===================== 主循环 ===============================
fprintf('\n开始 CPO-DBO（Permutation-safe）...\n');

for iter = 1:maxIter
   % t = iter / maxIter;
    t = min(1, iter / Tscale); 

    % 动态概率：前期开发+探索，后期增加逃跑（跳出局部）
    defenseProb = defenseProb0 * (1 - 0.6*t);
    attackProb  = attackProb0  * (1 - 0.2*t);
    escapeProb  = min(0.60, escapeProb0 + 0.55*t);

    % 归一化，确保三者和为1
    sProb = defenseProb + attackProb + escapeProb;
    defenseProb = defenseProb / sProb;
    attackProb  = attackProb  / sProb;
    escapeProb  = escapeProb  / sProb;

    % 精英保留：把当前最优的若干个体直接带入下一代（稳定性显著提升）
    [~, idxSort] = sort(fitness, 'ascend');
    eliteIdx = idxSort(1:eliteKeep);
    newPopulation = population;
    newFitness = fitness;
    for e = 1:eliteKeep
        newPopulation{e} = population{eliteIdx(e)};
        newFitness(e) = fitness(eliteIdx(e));
    end

    % 依次更新剩余个体
    for i = eliteKeep+1:M
        currentPath = population{i};
        currentFit  = fitness(i);

        % 依据群体类型设定“倾向”，但最终由动态概率决定行为
        if i <= popBall
            groupBias = 0.10; % 滚球群更偏防御
        elseif i <= popBall + popBreed
            groupBias = -0.05; % 繁殖群更偏攻击
        else
            groupBias = -0.10; % 小蜣螂更偏逃跑
        end

        % 选择行为
        r = rand();
        pD = min(0.95, max(0.05, defenseProb + groupBias));
        pA = min(0.95, max(0.05, attackProb  - groupBias/2));
        % 重新归一化
        pS = pD + pA + escapeProb;
        pD = pD / pS;
        pA = pA / pS;
        pE = escapeProb / pS;

        if r < pD
            % ---------- 防御（开发：向gBest学习 + 轻微变异 + 2-opt） ----------
            newPath = defenseBehavior_perm(currentPath, gBestPath, numNodes);
            newPath = mutate_perm(newPath, numNodes, 0.60, 0.15, 0.15);  % swap/insert/inversion 概率
            newPath = twoOptLocal(newPath, doseMatrix, twoOptIters_defense);
        elseif r < pD + pA
            % ---------- 攻击（探索：与pBest/随机体 OX交叉 + 变异 + 2-opt） ----------
            mateIdx = randi(M);
            matePath = pBestPath{mateIdx};
            newPath = oxCrossover(currentPath, matePath);
            newPath = mutate_perm(newPath, numNodes, 0.35, 0.20, 0.30);
            newPath = twoOptLocal(newPath, doseMatrix, twoOptIters_attack);
        else
            % ---------- 逃跑（跳出局部：较强扰动 + (可选)少量2-opt） ----------
            newPath = escapeBehavior_perm(currentPath, numNodes);
            newPath = mutate_perm(newPath, numNodes, 0.25, 0.25, 0.45);
            newPath = twoOptLocal(newPath, doseMatrix, twoOptIters_escape);
        end

        % 确保可行（这里应当几乎不触发修复，但保底）
        newPath = validatePath_det(newPath, numNodes);

        newFit = calculateFitness(newPath, doseMatrix);

        % 选择：用贪婪更新（也可以用 Metropolis，但先保证稳定）
        if newFit < currentFit
            newPopulation{i} = newPath;
            newFitness(i) = newFit;

            % 更新个体最优
            if newFit < pBestFitness(i)
                pBestPath{i} = newPath;
                pBestFitness(i) = newFit;
            end
        else
            % 不更新则保留原个体（稳定）
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
        fprintf('Iter %d/%d | Best dose: %.4f μSv | (pD=%.2f, pA=%.2f, pE=%.2f)\n',...
            iter, maxIter, gBestFitness, defenseProb, attackProb, escapeProb);
    end
end

%% ===================== 可视化 ===============================
% 收敛曲线
figure('Name', '算法收敛曲线', 'Position', [100, 100, 800, 400]);
plot(convergence, 'LineWidth', 2);
xlabel('迭代次数');
ylabel('最优剂量 (μSv)');
title('CPO-DBO（Permutation-safe）收敛曲线 (案例3: 51节点)');
grid on;

% 节点与放射源分布
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
    text(xCoords(i)+offset, yCoords(i)+offset/2, sprintf('%d', nodeIDs(i)), ...
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

% 最优路径图
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
    text(xCoords(i)+offset, yCoords(i)+offset/2, sprintf('%d', nodeIDs(i)), ...
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
    if mod(k, 6) == 0
        midX = (xCoords(a) + xCoords(b)) / 2;
        midY = (yCoords(a) + yCoords(b)) / 2;
        text(midX, midY, sprintf('%d→%d', a, b), ...
            'FontSize', 8, 'BackgroundColor', [1,1,1,0.7], ...
            'HorizontalAlignment', 'center');
    end
end

plot(xCoords(1), yCoords(1), 'o', 'MarkerSize', 12, ...
    'MarkerFaceColor', 'g', 'MarkerEdgeColor', 'k');
text(xCoords(1), yCoords(1)+5, '起点', 'FontSize', 10, ...
    'FontWeight', 'bold', 'Color', 'g', 'HorizontalAlignment', 'center');

legend('节点位置', '放射源位置', '路径', 'Location', 'northeast', 'FontSize', 12);
box on; grid minor;
set(gca, 'GridColor', [0.85, 0.85, 0.85], 'GridAlpha', 0.7);
hold off;

% 节点间剂量热力图
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

doseValues = doseMatrix;
doseValues(doseValues == 0) = NaN;
p95 = prctile(allDoses, 95);
for i = 1:numNodes
    for j = 1:numNodes
        if ~isnan(doseValues(i,j)) && doseValues(i,j) > p95
            text(j, i, sprintf('%.0f', doseValues(i,j)), ...
                 'HorizontalAlignment', 'center', 'FontSize', 6, ...
                 'FontWeight', 'bold', 'Color', 'k');
        end
    end
end

for i = 1:size(hotSpots, 1)
    nodeA = hotSpots(i, 1);
    nodeB = hotSpots(i, 2);
    rectangle('Position', [nodeB-0.5, nodeA-0.5, 1, 1], ...
              'EdgeColor', 'r', 'LineWidth', 1.5, 'LineStyle', '--');
end

annotation('textbox', [0.15, 0.02, 0.7, 0.05], ...
           'String', '注: 红色虚线框表示高辐射边（前5%），数值仅显示 >95% 分位数的剂量', ...
           'FontSize', 11, 'EdgeColor', 'none', 'HorizontalAlignment', 'center');

axis square;
set(gcf, 'Color', 'w');
hold off;

%% 输出结果
fprintf('\n===== 案例3复现结果 (Permutation-safe CPO-DBO) =====\n');
fprintf('最优路径: %s\n', strjoin(cellstr(num2str(gBestPath')), '->'));
fprintf('最优路径总剂量: %.4f μSv\n', gBestFitness);
fprintf('算法总迭代次数: %d\n', maxIter);

fprintf('放射源位置:\n');
for i = 1:size(sources, 1)
    fprintf('放射源%d: (%.1f, %.1f)\n', i, sources(i,1), sources(i,2));
end

%% ===================== 函数区 ==============================

% ---------- 防御行为（Permutation-safe）：向gBest学习 + 小扰动 ----------
function newPath = defenseBehavior_perm(path, bestPath, numNodes)
    % 取出排列部分（去掉首尾1）
    p = path(2:end-1);
    b = bestPath(2:end-1);

    % 方式1：对齐一小段（采用“位置学习”，但用 OX 保证可行）
    % 将b当作mate，用OX交叉能更自然地继承顺序信息
    if rand() < 0.80
        child = oxPerm(p, b);
    else
        child = p;
    end

    % 方式2：轻微扰动（swap）
    if rand() < 0.50
        child = swapPerm(child);
    end

    newPath = [1, child, 1];
    newPath = validatePath_det(newPath, numNodes);
end

% ---------- 逃跑行为（Permutation-safe）：较强扰动 ----------
function newPath = escapeBehavior_perm(path, numNodes)
    p = path(2:end-1);

    % 强扰动：随机打乱一段（scramble）+ inversion
    if rand() < 0.80
        p = scramblePerm(p);
    end
    if rand() < 0.60
        p = inversionPerm(p);
    end

    newPath = [1, p, 1];
    newPath = validatePath_det(newPath, numNodes);
end

% ---------- OX交叉（Path版本，含首尾1） ----------
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

    % 随机选择片段
    i = randi([1, n-1]);
    j = randi([i+1, n]);

    child(i:j) = a(i:j);

    % 按b的顺序填充剩余位置
    fillPos = [1:i-1, j+1:n];
    bSeq = b(~ismember(b, child(i:j)));

    child(fillPos) = bSeq;
end

% ---------- 变异组合：swap / insert / inversion ----------
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

% ---------- 2-opt 局部改进（少量迭代即可显著稳定） ----------
function path2 = twoOptLocal(path, doseMatrix, iters)
    if iters <= 0
        path2 = path;
        return;
    end
    p = path;
    n = numel(p);
    % p 结构：[1 ... 1]，中间长度为 numNodes-1
    for k = 1:iters
        i = randi([2, n-3]);   % 避开首尾1
        j = randi([i+1, n-2]);
        if j <= i
            continue;
        end

        old = edgeCost(p(i-1), p(i), doseMatrix) + edgeCost(p(j), p(j+1), doseMatrix);
        new = edgeCost(p(i-1), p(j), doseMatrix) + edgeCost(p(i), p(j+1), doseMatrix);

        if new < old
            p(i:j) = fliplr(p(i:j));
        end
    end
    path2 = p;
end

function c = edgeCost(a, b, doseMatrix)
    c = doseMatrix(a,b);
end

% ---------- 适应度：路径总剂量 ----------
function totalDose = calculateFitness(path, doseMatrix)
    totalDose = 0;
    for i = 1:(length(path)-1)
        totalDose = totalDose + doseMatrix(path(i), path(i+1));
    end
end

% ---------- 可行性校验（确定性、不会随机重置） ----------
function validPath = validatePath_det(path, numNodes)
    % 目标：确保 path = [1, perm(2..numNodes), 1] 且 perm无重复无缺失

    % 强制首尾为1
    path(1) = 1;
    path(end) = 1;

    mid = path(2:end-1);
    % 去掉非法（<2 或 >numNodes 或 =1）
    mid = mid(mid>=2 & mid<=numNodes);

    % 去重（保留第一次出现）
    [~, ia] = unique(mid, 'stable');
    mid = mid(sort(ia));

    % 补齐缺失节点（按节点编号顺序补齐，确定性）
    allNodes = 2:numNodes;
    missing = setdiff(allNodes, mid, 'stable');
    mid = [mid, missing];

    % 截断到正确长度
    mid = mid(1:(numNodes-1));

    validPath = [1, mid, 1];
end