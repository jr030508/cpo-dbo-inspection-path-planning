%% 清理环境
clc;
clear;
close all;

%% 节点坐标数据（表4）
% 论文中表4: 30个节点的位置（单位：米）
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

%% 假设放射源位置（基于图5）
% 根据论文图5，假设4个放射源位于高辐射区域
sources = [
    20, 70;  % 放射源1 (左上区域)
    40, 40;  % 放射源2 (中心区域)
    60, 60;  % 放射源3 (右上区域)
    70, 20;  % 放射源4 (右下区域)
];

%% 计算剂量矩阵
% 基于节点到放射源的距离计算剂量
doseMatrix = zeros(numNodes, numNodes);

% 计算每个节点的辐射水平
nodeRadLevel = zeros(numNodes, 1);
for i = 1:numNodes
    totalRad = 0;
    for s = 1:size(sources, 1)
        dist = norm([xCoords(i)-sources(s,1), yCoords(i)-sources(s,2)]);
        totalRad = totalRad + 1000 / (dist^2 + 5); % 辐射强度随距离衰减
    end
    nodeRadLevel(i) = totalRad;
end

% 计算节点间路径的剂量
for i = 1:numNodes
    for j = i+1:numNodes
        % 计算路径距离
        dist = norm([xCoords(i)-xCoords(j), yCoords(i)-yCoords(j)]);
        
        % 计算路径上的平均辐射水平
        avgRad = (nodeRadLevel(i) + nodeRadLevel(j)) / 2;
        
        % 考虑路径上可能经过高辐射区域
        midX = (xCoords(i) + xCoords(j)) / 2;
        midY = (yCoords(i) + yCoords(j)) / 2;
        midRad = 0;
        for s = 1:size(sources, 1)
            distToSource = norm([midX-sources(s,1), midY-sources(s,2)]);
            midRad = midRad + 500 / (distToSource^2 + 5);
        end
        
        % 最终剂量 = 距离 × (平均辐射水平 + 中点辐射水平)
        dose = dist * (avgRad + midRad) / 100;
        
        doseMatrix(i,j) = dose;
        doseMatrix(j,i) = dose;
    end
end

% 添加高辐射热点（基于图5）
hotSpots = [6, 20; 7, 20; 20, 26; 26, 21]; % 高辐射节点对
for i = 1:size(hotSpots, 1)
    nodeA = hotSpots(i, 1);
    nodeB = hotSpots(i, 2);
    doseMatrix(nodeA, nodeB) = doseMatrix(nodeA, nodeB) * 5; % 显著增加剂量
    doseMatrix(nodeB, nodeA) = doseMatrix(nodeB, nodeA) * 5;
end

% 对角线设为零
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

% 将剂量矩阵保存到CSV文件
csvwrite('dose_matrix_case2.csv', doseMatrix);
fprintf('\n暴露代价矩阵已保存到文件: dose_matrix_case2.csv\n');

%% CPO算法参数设置
M = 200;              % 冠豪猪种群数量
maxIter = 100;       % 最大迭代次数

% CPO特定参数
defenseProb = 0.7;    % 防御行为概率
attackProb = 0.3;     % 攻击行为概率
escapeProb = 0.2;     % 逃跑行为概率
quillEffect = 0.7;    % 刺毛效应强度

%% 初始化冠豪猪种群
% 每个冠豪猪代表一个路径排列
population = cell(M, 1);
fitness = zeros(M, 1);

% 初始化路径
for i = 1:M
    path = [1, randperm(numNodes-1)+1, 1]; % 起点和终点为1
    population{i} = path;
    fitness(i) = calculateFitness(path, doseMatrix);
end

% 初始化个体最优和全局最优
pBestPath = population;
pBestFitness = fitness;
[gBestFitness, gBestIdx] = min(fitness);
gBestPath = population{gBestIdx};

% 记录收敛曲线
convergence = zeros(maxIter, 1);
convergence(1) = gBestFitness;

%% CPO优化主循环
fprintf('\n开始冠豪猪优化(CPO)算法...\n');
for iter = 1:maxIter
    % 动态参数更新
    currentDefenseProb = defenseProb * (1 - iter/maxIter);  % 防御概率随迭代减少
    currentAttackProb = attackProb * (1 - iter/maxIter);    % 攻击概率随迭代减少
    currentEscapeProb = escapeProb + (0.3 * iter/maxIter);  % 逃跑概率随迭代增加
    
    for i = 1:M
        currentPath = population{i};
        currentFit = fitness(i);
        
        % ===== 行为选择 =====
        randBehavior = rand();
        
        if randBehavior < currentDefenseProb
            %% 防御行为（刺毛防御） - 局部优化
            newPath = defenseBehavior(currentPath, gBestPath, numNodes, quillEffect);
            
        elseif randBehavior < currentDefenseProb + currentAttackProb
            %% 攻击行为（刺毛攻击） - 探索新区域
            newPath = attackBehavior(currentPath, pBestPath, M, numNodes);
            
        else
            %% 逃跑行为 - 大幅扰动跳出局部最优
            newPath = escapeBehavior(currentPath, numNodes);
        end
        
        % 验证新路径
        newPath = validatePath(newPath, numNodes);
        newFit = calculateFitness(newPath, doseMatrix);
        
        % 更新当前个体
        if newFit < currentFit
            population{i} = newPath;
            fitness(i) = newFit;
            
            % 更新个体最优
            if newFit < pBestFitness(i)
                pBestPath{i} = newPath;
                pBestFitness(i) = newFit;
            end
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
    if mod(iter, 50) == 0
        fprintf('CPO迭代 %d/%d, 当前最优暴露代价: %.2f a.u.\n', iter, maxIter, gBestFitness);
    end
end

%% 结果可视化
% 绘制收敛曲线
figure('Name', '算法收敛曲线', 'Position', [100, 100, 800, 400]);
plot(convergence, 'LineWidth', 2);
xlabel('迭代次数');
ylabel('最优暴露代价 (a.u.)');
title('冠豪猪优化(CPO)收敛曲线 (案例2)');
grid on;
legend('最优暴露代价', 'Location', 'northeast');

%% ===== 增强节点图 =====
figure('Name', '节点与放射源分布', 'Position', [100, 100, 900, 700]);
hold on;

% 创建背景网格
grid on;
axis equal;
xlim([0 80]);
ylim([0 80]);
set(gca, 'FontSize', 12);
xlabel('X坐标 (m)', 'FontSize', 14);
ylabel('Y坐标 (m)', 'FontSize', 14);
title('节点与放射源分布图', 'FontSize', 16);

% 绘制节点 - 使用不同颜色和大小
scatter(xCoords, yCoords, 200, 'filled', 'MarkerFaceColor', [0.1, 0.5, 0.8], ...
    'MarkerEdgeColor', [0, 0.2, 0.4], 'LineWidth', 1.5);

% 添加清晰的节点标签
offset = 1.5; % 文本偏移量
for i = 1:numNodes
    % 动态调整文本位置避免重叠
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

% 标记放射源 - 使用醒目标记
scatter(sources(:,1), sources(:,2), 300, 'r', 'pentagram', 'filled', 'LineWidth', 1.5);
text(sources(:,1), sources(:,2)-3, {'放射源1','放射源2','放射源3','放射源4'},...
    'FontSize', 14, 'FontWeight', 'bold', 'Color', 'r', ...
    'HorizontalAlignment', 'center', 'VerticalAlignment', 'top');

% 添加图例
legend('节点位置', '放射源位置', 'Location', 'northeast', 'FontSize', 12);

% 添加网格和边框
box on;
grid minor;
set(gca, 'GridColor', [0.85, 0.85, 0.85], 'GridAlpha', 0.7);

hold off;

%% ===== 最优路径图 =====
figure('Name', '最优路径规划', 'Position', [100, 100, 900, 700]);
hold on;

% 创建背景
grid on;
axis equal;
xlim([0 80]);
ylim([0 80]);
set(gca, 'FontSize', 12);
xlabel('X坐标 (m)', 'FontSize', 14);
ylabel('Y坐标 (m)', 'FontSize', 14);
title(sprintf('案例2最优路径规划 (总暴露代价: %.2f a.u.)', gBestFitness), 'FontSize', 16);

% 绘制节点
scatter(xCoords, yCoords, 200, 'filled', 'MarkerFaceColor', [0.1, 0.5, 0.8], ...
    'MarkerEdgeColor', [0, 0.2, 0.4], 'LineWidth', 1.5);

% 添加节点标签
for i = 1:numNodes
    % 动态调整文本位置避免重叠
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

% 标记放射源
scatter(sources(:,1), sources(:,2), 300, 'r', 'pentagram', 'filled', 'LineWidth', 1.5);
text(sources(:,1), sources(:,2)-3, {'放射源1','放射源2','放射源3','放射源4'},...
    'FontSize', 14, 'FontWeight', 'bold', 'Color', 'r', ...
    'HorizontalAlignment', 'center', 'VerticalAlignment', 'top');

% 绘制最优路径
pathColors = jet(length(gBestPath)-1); % 使用不同颜色表示路径段
for i = 1:(length(gBestPath)-1)
    startIdx = gBestPath(i);
    endIdx = gBestPath(i+1);
    
    % 绘制路径线
    plot([xCoords(startIdx), xCoords(endIdx)],...
         [yCoords(startIdx), yCoords(endIdx)],...
         'Color', pathColors(i,:), 'LineWidth', 3, 'LineStyle', '-');
    
    % 添加路径方向标记
    midX = (xCoords(startIdx) + xCoords(endIdx)) / 2;
    midY = (yCoords(startIdx) + yCoords(endIdx)) / 2;
    text(midX, midY, sprintf('%d→%d', startIdx, endIdx), ...
        'FontSize', 10, 'BackgroundColor', [1,1,1,0.7], ...
        'HorizontalAlignment', 'center');
end

% 添加起点和终点标记
plot(xCoords(1), yCoords(1), 'o', 'MarkerSize', 15, 'MarkerFaceColor', 'g', 'MarkerEdgeColor', 'k');
text(xCoords(1), yCoords(1)+3, '起点', 'FontSize', 12, 'FontWeight', 'bold', 'Color', 'g', ...
    'HorizontalAlignment', 'center');

% 添加图例
legend('节点位置', '放射源位置', '路径', 'Location', 'northeast', 'FontSize', 12);

% 添加网格和边框
box on;
grid minor;
set(gca, 'GridColor', [0.85, 0.85, 0.85], 'GridAlpha', 0.7);

hold off;

%% ===== 增强型剂量热力图 =====
figure('Name', '节点间暴露代价热力图', 'Position', [100, 100, 1000, 850]);

% 创建自定义颜色映射
colors = [0.0, 0.0, 1.0;   % 蓝色 (低辐射)
          0.0, 1.0, 0.0;   % 绿色
          1.0, 1.0, 0.0;   % 黄色
          1.0, 0.5, 0.0;   % 橙色
          1.0, 0.0, 0.0];  % 红色 (高辐射)
colorPos = [0, 0.3, 0.6, 0.8, 1];
customMap = interp1(colorPos, colors, linspace(0,1,256));

% 创建热力图
h = imagesc(doseMatrix);
colormap(customMap);
c = colorbar;
c.Label.String = '暴露代价 (a.u.)';
c.Label.FontSize = 12;
c.FontSize = 11;

% 设置坐标轴
set(gca, 'XTick', 1:numNodes, 'YTick', 1:numNodes, ...
         'FontSize', 10, 'XAxisLocation', 'top');
xlabel('目标节点', 'FontSize', 14, 'FontWeight', 'bold');
ylabel('起始节点', 'FontSize', 14, 'FontWeight', 'bold');
title('节点间暴露代价矩阵热力图 (案例2)', 'FontSize', 16, 'FontWeight', 'bold');

% 添加网格线
hold on;
for i = 1:numNodes+1
    plot([0.5, numNodes+0.5], [i-0.5, i-0.5], 'k-', 'LineWidth', 0.5);
    plot([i-0.5, i-0.5], [0.5, numNodes+0.5], 'k-', 'LineWidth', 0.5);
end

% 添加剂量值文本 - 只显示较高剂量的值
doseValues = doseMatrix;
doseValues(doseValues == 0) = NaN; % 不显示对角线

for i = 1:numNodes
    for j = 1:numNodes
        if doseValues(i,j) > 100
            text(j, i, sprintf('%.0f', doseValues(i,j)), ...
                 'HorizontalAlignment', 'center', ...
                 'FontSize', 9, ...
                 'FontWeight', 'bold', ...
                 'Color', 'k');
        end
    end
end

% 添加高辐射标记
for i = 1:size(hotSpots, 1)
    nodeA = hotSpots(i, 1);
    nodeB = hotSpots(i, 2);
    rectangle('Position', [nodeB-0.5, nodeA-0.5, 1, 1], ...
              'EdgeColor', 'r', 'LineWidth', 2, 'LineStyle', '--');
end

% 添加图例说明
annotation('textbox', [0.15, 0.02, 0.7, 0.05], ...
           'String', '注: 红色虚线框表示高辐射路径，数值仅显示大于100a.u.的剂量', ...
           'FontSize', 11, 'EdgeColor', 'none', 'HorizontalAlignment', 'center');

% 优化布局
axis square;
set(gcf, 'Color', 'w');
hold off;

%% ===== 辐射环境热力图 =====
figure('Name', '辐射环境热力图', 'Position', [100, 100, 900, 700]);
hold on;

% 创建网格
[xGrid, yGrid] = meshgrid(0:1:80, 0:1:80);
radGrid = zeros(size(xGrid));

% 计算网格点的辐射水平
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

% 绘制热力图
contourf(xGrid, yGrid, radGrid, 30, 'LineStyle', 'none');
colormap(jet);
c = colorbar;
c.Label.String = '辐射水平';
c.Label.FontSize = 12;
title('辐射环境热力图', 'FontSize', 16);
xlabel('X坐标 (m)', 'FontSize', 14);
ylabel('Y坐标 (m)', 'FontSize', 14);

% 绘制节点
scatter(xCoords, yCoords, 100, 'k', 'filled', 'MarkerEdgeColor', 'w', 'LineWidth', 1.5);

% 添加节点标签
for i = 1:numNodes
    text(xCoords(i)+1.5, yCoords(i)+1.5, sprintf('%d', nodeIDs(i)), ...
         'FontSize', 11, 'FontWeight', 'bold', 'Color', 'w', ...
         'BackgroundColor', [0,0,0,0.5], ...
         'HorizontalAlignment', 'center');
end

% 标记放射源
scatter(sources(:,1), sources(:,2), 200, 'r', 'pentagram', 'filled', 'LineWidth', 1.5);

% 添加等高线
contour(xGrid, yGrid, radGrid, 10, 'LineColor', 'k', 'LineWidth', 0.5);

% 添加图例
legend('节点位置', '放射源位置', 'Location', 'northeast', 'FontSize', 12);

% 设置坐标轴
axis equal;
xlim([0 80]);
ylim([0 80]);

hold off;

%% 输出结果
fprintf('\n===== 案例2复现结果 (CPO算法) =====\n');
fprintf('最优路径: %s\n', strjoin(cellstr(num2str(gBestPath')), '->'));
fprintf('最优路径总暴露代价: %.2f a.u.\n', gBestFitness);
fprintf('算法总迭代次数: %d\n', maxIter);

%% CPO行为函数
function newPath = defenseBehavior(path, bestPath, numNodes, quillEffect)
    % 刺毛防御行为 - 局部优化
    newPath = path;
    
    % 向全局最优学习
    if rand() < quillEffect
        % 随机选择一段最优路径片段
        segLength = randi([2, min(5, numNodes-2)]);
        startIdx = randi([2, numNodes-segLength]);
        bestSeg = bestPath(startIdx:startIdx+segLength-1);
        
        % 替换当前路径的对应片段
        newPath(startIdx:startIdx+segLength-1) = bestSeg;
    end
    
    % 局部扰动
    if rand() < 0.6
        % 随机交换两个节点
        swapIdx = randperm(numNodes-2, 2) + 1; % 避开起点终点
        temp = newPath(swapIdx(1));
        newPath(swapIdx(1)) = newPath(swapIdx(2));
        newPath(swapIdx(2)) = temp;
    end
end

function newPath = attackBehavior(path, bestPaths, popSize, numNodes)
    % 刺毛攻击行为 - 探索新区域
    newPath = path;
    
    % 随机选择另一个个体
    otherIdx = randi(popSize);
    otherPath = bestPaths{otherIdx};
    
    % 交叉操作
    if rand() < 0.7
        % 随机选择交叉点
        crossPoint = randi([2, numNodes-2]);
        
        % 替换后半部分
        newPath(crossPoint:end-1) = otherPath(crossPoint:end-1);
    end
    
    % 大幅变异
    if rand() < 0.4
        % 随机反转一段路径
        segStart = randi([2, numNodes-3]);
        segEnd = segStart + randi([1, min(4, numNodes-segStart-1)]);
        newPath(segStart:segEnd) = fliplr(newPath(segStart:segEnd));
    end
end

function newPath = escapeBehavior(path, numNodes)
    % 逃跑行为 - 大幅扰动跳出局部最优
    newPath = path;
    
    % 随机重新排序中间节点
    middleNodes = path(2:end-1);
    newPath(2:end-1) = middleNodes(randperm(length(middleNodes)));
    
    % 随机插入新节点（确保所有节点都被访问）
    if rand() < 0.3
        missingNodes = setdiff(1:numNodes, newPath);
        if ~isempty(missingNodes)
            insertPos = randi([2, length(newPath)-1]);
            newPath = [newPath(1:insertPos-1), missingNodes(1), newPath(insertPos:end)];
        end
    end
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
    [~, uniqueIdx] = unique(path(2:end-1), 'stable');
    validPath = [1, path(uniqueIdx+1), 1];
    
    % 确保长度正确
    if length(validPath) ~= numNodes + 1
        validPath = [1, randperm(numNodes-1) + 1, 1];
    end
end