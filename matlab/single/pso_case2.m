%% 核设施辐射环境路径规划 - 案例1完整复现（最终版）
% 作者: Zhuang Wang, Jiejin Cai (论文)
% MATLAB实现: 基于改进粒子群算法(混沌PSO)
% 日期: 2023-10-15

%% 清理环境
clc;
clear;
close all;

%% 节点坐标数据（表1）
nodeData = [
    1, 30, 0;
    2, 30, 10;
    3, 30, 20;
    4, 35, 27;
    5, 26, 17;
    6, 16, 26;
    7, 21, 21;
    8, 10, 21;
    9, 2, 25;
    10, 1, 8;
    11, 10, 14;
    12, 21, 13;
    13, 25, 24;
    14, 3, 14;
];

nodeIDs = nodeData(:,1);
xCoords = nodeData(:,2);
yCoords = nodeData(:,3);
numNodes = length(nodeIDs);

%% 有效剂量数据（表2）
% 论文中表2: 14个节点间的暴露代价(a.u.)
doseMatrix = [
999, 4.88, 13.81, 39.27, 9.24, 327.36, 29.92, 28.4, 682.79, 6.08, 13.62, 6.83, 13.97, 10.04;
4.88, 999, 8.95, 21.03, 6.22, 24.66, 17.79, 35.36, 37.13, 11.77, 21.83, 7, 10.5, 20.61;
13.81, 8.95, 999, 4.85, 3.84, 7.88, 7.39, 18.06, 15.05, 79.24, 179.93, 12.67, 3.51, 1362.88;
39.27, 21.03, 4.85, 999, 7.91, 5.82, 7.56, 11.97, 9.17, 51.93, 66.4, 17.89, 3.67, 39.98;
9.24, 6.22, 3.84, 7.91, 999, 11.95, 8.81, 39.9, 32.32, 34.1, 66.56, 8.96, 4.97, 117.13;
327.36, 24.66, 7.88, 5.82, 11.95, 999, 4.04, 3.45, 3.26, 13.89, 45.57, 29.75, 3.55, 8.97;
29.92, 17.79, 7.39, 7.56, 8.81, 4.04, 999, 9.07, 9.14, 68.66, 73.62, 3415.27, 3.38, 36.07;
28.4, 35.36, 18.06, 11.97, 39.9, 3.45, 9.07, 999, 2.59, 7.12, 14.88, 20.38, 8.94, 4.99;
682.79, 37.13, 15.05, 9.17, 32.32, 3.26, 9.14, 2.59, 999, 2.81, 8.94, 26.21, 7.71, 2.17;
6.08, 11.77, 79.24, 51.93, 34.1, 13.89, 68.66, 7.12, 2.81, 999, 4.94, 11.32, 88.61, 1.11;
13.62, 21.83, 179.93, 66.4, 66.56, 45.57, 73.62, 14.88, 8.94, 4.94, 999, 18.52, 83.8, 4.69;
6.83, 7, 12.67, 17.89, 8.96, 29.75, 3415.27, 20.38, 26.21, 11.32, 18.52, 999, 29.85, 20.55;
13.97, 10.5, 3.51, 3.67, 4.97, 3.55, 3.38, 8.94, 7.71, 88.61, 83.8, 29.85, 999, 28.9;
10.04, 20.61, 1362.88, 39.98, 117.13, 8.97, 36.07, 4.99, 2.17, 1.11, 4.69, 20.55, 28.9, 999;
];

% 将剂量矩阵对角线设为0（自身到自身）
%doseMatrix(logical(eye(numNodes))) = 0;

%% 算法参数设置（根据论文）
M = 200;             % 粒子数量
maxChaosIter = 20;  % 混沌优化迭代次数
maxPSOIter = 2000;   % PSO总迭代次数

%% 初始化粒子群
% 初始化混沌变量
chaosPop = initChaosPopulation(M, numNodes);

% 将混沌变量转换为路径
population = cell(M, 1);
fitness = zeros(M, 1);
for i = 1:M
    path = chaosToPath(chaosPop(i, :), numNodes);
    population{i} = path;
    fitness(i) = calculateFitness(path, doseMatrix);
end

% 初始化个体最优和全局最优
pBestPath = population;
pBestFitness = fitness;
[gBestFitness, gBestIdx] = min(fitness);
gBestPath = population{gBestIdx};

% 记录收敛曲线
convergence = zeros(maxChaosIter + maxPSOIter, 1);
convergence(1) = gBestFitness;

%% 混沌优化阶段（论文步骤1-3）
fprintf('开始混沌优化阶段...\n');
for iter = 1:maxChaosIter
    % 更新混沌变量
    chaosPop = updateChaos(chaosPop);
    
    % 更新粒子位置和适应度
    for i = 1:M
        newPath = chaosToPath(chaosPop(i, :), numNodes);
        newFitness = calculateFitness(newPath, doseMatrix);
        
        % 更新个体最优
        if newFitness < pBestFitness(i)
            pBestPath{i} = newPath;
            pBestFitness(i) = newFitness;
        end
    end
    
    % 更新全局最优
    [minFitness, minIdx] = min(pBestFitness);
    if minFitness < gBestFitness
        gBestFitness = minFitness;
        gBestPath = pBestPath{minIdx};
    end
    
    % 记录收敛曲线
    convergence(iter) = gBestFitness;
    
    % 显示迭代信息
    fprintf('混沌迭代 %d/%d, 当前最优暴露代价: %.2f a.u.\n', iter, maxChaosIter, gBestFitness);
end

%% PSO优化阶段（论文步骤4-5）
fprintf('\n开始PSO优化阶段...\n');
for iter = 1:maxPSOIter
    % 对每个粒子进行更新
    for i = 1:M
        % 获取当前粒子的路径
        currentPath = population{i};
        
        % 应用交叉操作（F1和F2）
        crossPath = crossoverOperator(currentPath, pBestPath{i}, gBestPath, numNodes);
        
        % 应用变异操作（F3）
        newPath = mutationOperator(crossPath);
        
        % 确保新路径有效
        newPath = validatePath(newPath, numNodes);
        newFitness = calculateFitness(newPath, doseMatrix);
        
        % 更新个体最优
        if newFitness < pBestFitness(i)
            pBestPath{i} = newPath;
            pBestFitness(i) = newFitness;
        end
        
        % 更新当前粒子
        population{i} = newPath;
        fitness(i) = newFitness;
    end
    
    % 更新全局最优
    [minFitness, minIdx] = min(pBestFitness);
    if minFitness < gBestFitness
        gBestFitness = minFitness;
        gBestPath = pBestPath{minIdx};
    end
    
    % 记录收敛曲线
    convergence(maxChaosIter + iter) = gBestFitness;
    
    % 显示迭代信息
    if mod(iter, 50) == 0
        fprintf('PSO迭代 %d/%d, 当前最优暴露代价: %.2f a.u.\n', iter, maxPSOIter, gBestFitness);
    end
end

%% 结果可视化
% 绘制收敛曲线
figure('Name', '算法收敛曲线', 'Position', [100, 100, 800, 400]);
plot(convergence, 'LineWidth', 2);
xlabel('迭代次数');
ylabel('最优暴露代价 (a.u.)');
title('混沌PSO收敛曲线');
grid on;
xline(maxChaosIter, 'r--', 'LineWidth', 1.5, 'Label', '混沌优化结束');
legend('最优暴露代价', 'Location', 'northeast');

% 绘制最优路径
figure('Name', '最优路径规划', 'Position', [100, 100, 800, 600]);
hold on;

% 绘制节点
scatter(xCoords, yCoords, 100, 'filled', 'MarkerFaceColor', [0.2,0.6,0.9]);
text(xCoords, yCoords, cellstr(num2str(nodeIDs)),...
    'VerticalAlignment','bottom','HorizontalAlignment','right','FontSize',12);

% 标记高辐射区域
highRadNodes = find(sum(doseMatrix > 100, 2) > 5); % 找出高辐射节点
scatter(xCoords(highRadNodes), yCoords(highRadNodes), 120,...
    'Marker','pentagram','MarkerFaceColor','r','MarkerEdgeColor','k');

% 绘制最优路径
for i = 1:(length(gBestPath)-1)
    startIdx = gBestPath(i);
    endIdx = gBestPath(i+1);
    plot([xCoords(startIdx), xCoords(endIdx)],...
         [yCoords(startIdx), yCoords(endIdx)],...
         'b-o','LineWidth',2,'MarkerSize',8);
end

% 添加标签和标题
title(sprintf('最优路径规划 (总暴露代价: %.2f a.u.)', gBestFitness));
xlabel('X坐标 (m)');
ylabel('Y坐标 (m)');
legend('普通节点', '高辐射节点', '路径', 'Location', 'best');
grid on;
axis equal;
hold off;

% 绘制剂量热力图
figure('Name','节点间暴露代价热力图','Position',[100,100,800,600]);
imagesc(doseMatrix);
colormap(jet);
colorbar;
title('节点间暴露代价矩阵热力图 (a.u.)');
xlabel('目标节点');
ylabel('起始节点');
set(gca,'XTick',1:numNodes,'YTick',1:numNodes);
axis square;

%% 输出结果
fprintf('\n===== 案例1复现结果 =====\n');
fprintf('最优路径: %s\n', strjoin(cellstr(num2str(gBestPath')), '->'));
fprintf('最优路径总暴露代价: %.2f a.u.\n', gBestFitness);
fprintf('算法总迭代次数: %d (混沌优化:%d + PSO:%d)\n',...
    maxChaosIter+maxPSOIter, maxChaosIter, maxPSOIter);

%% 辅助函数：初始化混沌种群
function chaosPop = initChaosPopulation(M, n)
    chaosPop = zeros(M, n);
    for i = 1:M
        % 生成(0,1)之间的随机数，排除特殊值
        temp = rand(1, n);
        while any(ismember(temp, [0, 0.25, 0.5, 0.75, 1]))
            temp = rand(1, n);
        end
        chaosPop(i, :) = temp;
    end
end

%% 辅助函数：混沌变量转路径
function path = chaosToPath(chaosVector, numNodes)
    [~, idx] = sort(chaosVector); % 获取排序索引
    path = idx; % 初始路径
    
    % 确保节点1在路径中
    if ~ismember(1, path)
        path(1) = 1;
    end
    
    % 确保路径以节点1开始和结束
    path = [1, setdiff(path, 1, 'stable'), 1];
    
    % 确保路径包含所有节点
    if length(unique(path)) < numNodes
        missing = setdiff(1:numNodes, path);
        path = [path(1:end-1), missing, 1];
    end
end

%% 辅助函数：更新混沌变量
function newChaos = updateChaos(oldChaos)
    newChaos = 4 * oldChaos .* (1 - oldChaos);
end

%% 辅助函数：计算路径适应度
function totalDose = calculateFitness(path, doseMatrix)
    totalDose = 0;
    for i = 1:(length(path)-1)
        startNode = path(i);
        endNode = path(i+1);
        
        % 确保索引在有效范围内
        if startNode < 1 || startNode > size(doseMatrix,1) || ...
           endNode < 1 || endNode > size(doseMatrix,2)
            totalDose = totalDose + 1e6; % 惩罚无效路径
        else
            totalDose = totalDose + doseMatrix(startNode, endNode);
        end
    end
end

%% 辅助函数：交叉操作
function child = crossoverOperator(parent1, parent2, gBest, numNodes)
    % 确保路径长度一致
    n = length(parent1);
    
    % 选择两个交叉点（避开起点和终点）
    crossPoints = sort(randperm(n-2, 2) + 1);
    startCross = crossPoints(1);
    endCross = crossPoints(2);
    
    % 从父代1复制片段
    segment = parent1(startCross:endCross);
    
    % 从父代2复制剩余节点
    remaining = parent2;
    remaining(ismember(remaining, segment)) = []; % 移除已存在的节点
    
    % 构建子代路径
    child = zeros(1, n);
    child(1) = 1; % 起点固定为1
    child(n) = 1; % 终点固定为1
    
    % 填充中间部分
    insertPos = 2;
    for i = 1:length(remaining)
        if remaining(i) ~= 1 % 跳过起点/终点
            child(insertPos) = remaining(i);
            insertPos = insertPos + 1;
            if insertPos > n-1
                break;
            end
        end
    end
    
    % 插入交叉片段
    child(startCross:endCross) = segment;
    
    % 10%概率使用全局最优片段
    if rand() < 0.1
        gbSegment = gBest(startCross:endCross);
        child(startCross:endCross) = gbSegment;
    end
end

%% 辅助函数：变异操作
function mutated = mutationOperator(path)
    % 随机交换两个位置（避开起点和终点）
    n = length(path);
    if n > 3 % 确保有足够的节点进行变异
        idx = randperm(n-2, 2) + 1; % 选择中间节点
        mutated = path;
        temp = mutated(idx(1));
        mutated(idx(1)) = mutated(idx(2));
        mutated(idx(2)) = temp;
    else
        mutated = path;
    end
    
    % 10%概率反转子路径
    if rand() < 0.1 && n > 4
        points = sort(randperm(n-2, 2) + 1);
        mutated(points(1):points(2)) = fliplr(mutated(points(1):points(2)));
    end
end

%% 辅助函数：验证路径
function validPath = validatePath(path, numNodes)
    % 确保路径包含所有节点
    if length(unique(path)) < numNodes
        % 添加缺失的节点
        missingNodes = setdiff(1:numNodes, path);
        for i = 1:length(missingNodes)
            insertPos = randi([2, length(path)-1]); % 随机插入位置
            path = [path(1:insertPos-1), missingNodes(i), path(insertPos:end)];
        end
    end
    
    % 确保起点和终点为1
    path(1) = 1;
    path(end) = 1;
    
    % 移除重复节点（起点和终点除外）
    uniquePath = unique(path(2:end-1), 'stable');
    validPath = [1, uniquePath, 1];
    
    % 确保长度正确
    if length(validPath) ~= numNodes + 1
        validPath = [1, randperm(numNodes-1) + 1, 1];
    end
end