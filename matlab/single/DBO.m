%% 清理环境
clc;
clear;
close all;

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

%% 设置7个固定放射源位置（基于论文案例3的图11布局）
% 这些位置基于案例3的辐射环境热力图大致分布
sources = [
    35, 35;    % 放射源1 - 左下区域
    30, 110;   % 放射源2 - 右下区域
    88, 70;    % 放射源3 - 中心偏左
    88, 142;  % 放射源4 - 中心偏右
    138, 107;   % 放射源5 - 左上区域
    160, 50;  % 放射源6 - 右上区域
    168, 180;  % 放射源7 - 中上区域
];

% 设置每个放射源的强度系数（基于论文中的辐射水平）
sourceStrength = [8000, 12000, 10000, 15000, 9000, 13000, 11000];

%% 计算剂量矩阵
% 基于节点到放射源的距离计算剂量
doseMatrix = zeros(numNodes, numNodes);

% 计算每个节点的辐射水平
nodeRadLevel = zeros(numNodes, 1);
for i = 1:numNodes
    totalRad = 0;
    for s = 1:size(sources, 1)
        dist = norm([xCoords(i)-sources(s,1), yCoords(i)-sources(s,2)]);
        % 使用对应的放射源强度
        totalRad = totalRad + sourceStrength(s) / (dist^2 + 10); % 辐射强度随距离衰减
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
            % 使用对应的放射源强度的一半
            midRad = midRad + (sourceStrength(s)/2) / (distToSource^2 + 10);
        end
        
        % 最终剂量 = 距离 × (平均辐射水平 + 中点辐射水平)
        dose = dist * (avgRad + midRad) / 100;
        
        doseMatrix(i,j) = dose;
        doseMatrix(j,i) = dose;
    end
end

% 添加高辐射热点（基于剂量矩阵自动识别）
% 找出剂量最高的5%的路径作为高辐射热点
allDoses = doseMatrix(:);
allDoses(allDoses == 0) = []; % 移除对角线上的零值
threshold = prctile(allDoses, 95); % 取前5%的剂量值作为阈值

hotSpots = [];
for i = 1:numNodes
    for j = i+1:numNodes
        if doseMatrix(i,j) > threshold
            hotSpots = [hotSpots; i, j];
            % 显著增加高辐射路径的剂量
            doseMatrix(i,j) = doseMatrix(i,j) * 3;
            doseMatrix(j,i) = doseMatrix(j,i) * 3;
        end
    end
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


%% 算法参数设置（根据论文）
M = 400;              % 蜣螂种群数量（原粒子数量2000过大，改为200）
maxIter = 30000;       % 最大迭代次数（混沌优化200+PSO1000→总1200次）

%% DBO参数设置
popBall = round(M * 0.5);     % 滚球蜣螂数量
popBreed = round(M * 0.25);   % 繁殖蜣螂数量
popBaby = M - popBall - popBreed; % 小蜣螂数量

%% 初始化蜣螂种群
% 每个蜣螂代表一个路径排列
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

%% DBO优化主循环
fprintf('\n开始蜣螂优化(DBO)算法...\n');
for iter = 1:maxIter
    % 动态参数更新
    k = 1.0 - (iter/maxIter);    % 光线强度系数（线性递减）
    R = 1.0 - (iter/maxIter);    % 繁殖区域收缩系数
    
    % ===== 1. 滚球蜣螂更新 =====
    for i = 1:popBall
        % 随机选择alpha（-1或1）
        alpha = 1;
        if rand() > 0.5
            alpha = -1;
        end
        
        % 扰动因子
        b = rand();
        
        % 位置更新（模拟滚球行为）
        if iter > 1
            % 添加旋转扰动
            theta = 2 * pi * rand();
            rotation = [cos(theta), -sin(theta); sin(theta), cos(theta)];
            
            % 路径扰动：随机交换节点
            newPath = pBestPath{i};
            swapIdx = randperm(numNodes-2, 2) + 1; % 避开起点终点
            newPath(swapIdx) = newPath(swapIdx([2,1]));
            
            % 向全局最优学习
            if rand() < 0.3
                segStart = randi([2, numNodes-3]);
                segLength = randi([2, min(4, numNodes-segStart)]);
                bestSeg = gBestPath(segStart:segStart+segLength);
                newPath(segStart:segStart+segLength) = bestSeg;
            end
            
            % 更新路径并计算适应度
            newPath = validatePath(newPath, numNodes);
            newFitness = calculateFitness(newPath, doseMatrix);
            
            % 更新个体最优
            if newFitness < pBestFitness(i)
                pBestPath{i} = newPath;
                pBestFitness(i) = newFitness;
            end
        end
    end
    
    % ===== 2. 繁殖蜣螂更新 =====
    for i = popBall+1:popBall+popBreed
        % 安全区域繁殖（向全局最优靠近）
        newPath = pBestPath{i};
        
        % 路径交叉：与全局最优路径交换片段
        if rand() < 0.6
            crossPoint = randi([2, numNodes-2]); % 交叉点
            newPath(crossPoint:end-1) = gBestPath(crossPoint:end-1);
        end
        
        % 路径变异：随机交换两个节点
        if rand() < 0.4
            swapIdx = randperm(numNodes-2, 2) + 1; % 避开起点终点
            newPath(swapIdx) = newPath(swapIdx([2,1]));
        end
        
        % 更新路径并计算适应度
        newPath = validatePath(newPath, numNodes);
        newFitness = calculateFitness(newPath, doseMatrix);
        
        % 更新个体最优
        if newFitness < pBestFitness(i)
            pBestPath{i} = newPath;
            pBestFitness(i) = newFitness;
        end
    end
    
    % ===== 3. 小蜣螂更新 =====
    for i = popBall+popBreed+1:M
        newPath = pBestPath{i};
        
        % 向最优个体学习（觅食行为）
        if rand() < 0.7
            % 随机插入节点
            moveNode = randi([2, numNodes-1]);
            newPos = randi([2, numNodes-1]);
            
            % 移动节点到新位置
            newPath(newPath == moveNode) = [];
            newPath = [newPath(1:newPos-1), moveNode, newPath(newPos:end)];
        end
        
        % 随机扰动
        if rand() < 0.3
            swapIdx = randperm(numNodes-2, 2) + 1; % 避开起点终点
            newPath(swapIdx) = newPath(swapIdx([2,1]));
        end
        
        % 更新路径并计算适应度
        newPath = validatePath(newPath, numNodes);
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
    if mod(iter, 50) == 0
        fprintf('DBO迭代 %d/%d, 当前最优暴露代价: %.2f a.u.\n', iter, maxIter, gBestFitness);
    end
end

%% 结果可视化
% 绘制收敛曲线
figure('Name', '算法收敛曲线', 'Position', [100, 100, 800, 400]);
plot(convergence, 'LineWidth', 2);
xlabel('迭代次数');
ylabel('最优暴露代价 (a.u.)');
title('改进蜣螂优化器(DBO)收敛曲线 (案例2)');
grid on;
legend('最优暴露代价', 'Location', 'northeast');

%% ===== 增强节点图 =====
figure('Name', '节点与放射源分布', 'Position', [100, 100, 900, 700]);
hold on;

% 创建背景网格
grid on;
axis equal;
xlim([0 200]);
ylim([0 200]);
set(gca, 'FontSize', 12);
xlabel('X坐标 (m)', 'FontSize', 14);
ylabel('Y坐标 (m)', 'FontSize', 14);
title('节点与放射源分布图 (案例3: 51节点)', 'FontSize', 16);

% 绘制节点 - 使用不同颜色和大小
scatter(xCoords, yCoords, 100, 'filled', 'MarkerFaceColor', [0.1, 0.5, 0.8], ...
    'MarkerEdgeColor', [0, 0.2, 0.4], 'LineWidth', 1.5);

% 添加清晰的节点标签
offset = 3; % 文本偏移量增加以适应更大范围
for i = 1:numNodes
    % 动态调整文本位置避免重叠
    textPos = [xCoords(i) + offset, yCoords(i) + offset/2];
    
    text(textPos(1), textPos(2), sprintf('%d', nodeIDs(i)), ...
        'FontSize', 10, 'FontWeight', 'bold', 'Color', 'k', ...
        'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle');
end

% 标记放射源 - 使用醒目标记
scatter(sources(:,1), sources(:,2), 200, 'r', 'pentagram', 'filled', 'LineWidth', 1.5);
for s = 1:size(sources, 1)
    text(sources(s,1), sources(s,2)-5, sprintf('放射源%d', s),...
        'FontSize', 12, 'FontWeight', 'bold', 'Color', 'r', ...
        'HorizontalAlignment', 'center', 'VerticalAlignment', 'top');
end

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
xlim([0 200]);
ylim([0 200]);
set(gca, 'FontSize', 12);
xlabel('X坐标 (m)', 'FontSize', 14);
ylabel('Y坐标 (m)', 'FontSize', 14);
title(sprintf('案例3最优路径规划 (总暴露代价: %.2f a.u.)', gBestFitness), 'FontSize', 16);

% 绘制节点
scatter(xCoords, yCoords, 100, 'filled', 'MarkerFaceColor', [0.1, 0.5, 0.8], ...
    'MarkerEdgeColor', [0, 0.2, 0.4], 'LineWidth', 1.5);

% 添加节点标签
for i = 1:numNodes
    textPos = [xCoords(i) + offset, yCoords(i) + offset/2];
    
    text(textPos(1), textPos(2), sprintf('%d', nodeIDs(i)), ...
        'FontSize', 10, 'FontWeight', 'bold', 'Color', 'k', ...
        'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle');
end

% 标记放射源
scatter(sources(:,1), sources(:,2), 200, 'r', 'pentagram', 'filled', 'LineWidth', 1.5);
for s = 1:size(sources, 1)
    text(sources(s,1), sources(s,2)-5, sprintf('放射源%d', s),...
        'FontSize', 12, 'FontWeight', 'bold', 'Color', 'r', ...
        'HorizontalAlignment', 'center', 'VerticalAlignment', 'top');
end

% 绘制最优路径
pathColors = jet(length(gBestPath)-1); % 使用不同颜色表示路径段
for i = 1:(length(gBestPath)-1)
    startIdx = gBestPath(i);
    endIdx = gBestPath(i+1);
    
    % 绘制路径线
    plot([xCoords(startIdx), xCoords(endIdx)],...
         [yCoords(startIdx), yCoords(endIdx)],...
         'Color', pathColors(i,:), 'LineWidth', 2, 'LineStyle', '-');
    
    % 添加路径方向标记（仅显示部分，避免过于拥挤）
    if mod(i, 5) == 0 % 每5条路径显示一个标记
        midX = (xCoords(startIdx) + xCoords(endIdx)) / 2;
        midY = (yCoords(startIdx) + yCoords(endIdx)) / 2;
        text(midX, midY, sprintf('%d→%d', startIdx, endIdx), ...
            'FontSize', 8, 'BackgroundColor', [1,1,1,0.7], ...
            'HorizontalAlignment', 'center');
    end
end

% 添加起点和终点标记
plot(xCoords(1), yCoords(1), 'o', 'MarkerSize', 12, 'MarkerFaceColor', 'g', 'MarkerEdgeColor', 'k');
text(xCoords(1), yCoords(1)+5, '起点', 'FontSize', 10, 'FontWeight', 'bold', 'Color', 'g', ...
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
         'FontSize', 8, 'XAxisLocation', 'top');
xlabel('目标节点', 'FontSize', 14, 'FontWeight', 'bold');
ylabel('起始节点', 'FontSize', 14, 'FontWeight', 'bold');
title('节点间暴露代价矩阵热力图 (案例3: 51节点)', 'FontSize', 16, 'FontWeight', 'bold');

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
        if doseValues(i,j) > prctile(allDoses, 95) % 只显示前10%的高剂量值
            text(j, i, sprintf('%.0f', doseValues(i,j)), ...
                 'HorizontalAlignment', 'center', ...
                 'FontSize', 6, ...
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
              'EdgeColor', 'r', 'LineWidth', 1.5, 'LineStyle', '--');
end

% 添加图例说明
annotation('textbox', [0.15, 0.02, 0.7, 0.05], ...
           'String', '注: 红色虚线框表示高辐射路径，数值仅显示大于90%分位数的剂量', ...
           'FontSize', 11, 'EdgeColor', 'none', 'HorizontalAlignment', 'center');

% 优化布局
axis square;
set(gcf, 'Color', 'w');
hold off;

%% ===== 辐射环境热力图 =====
figure('Name', '辐射环境热力图', 'Position', [100, 100, 900, 700]);
hold on;

% 创建网格
[xGrid, yGrid] = meshgrid(0:2:200, 0:2:200); % 网格步长增加以适应更大范围
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
title('辐射环境热力图 (案例3: 51节点)', 'FontSize', 16);
xlabel('X坐标 (m)', 'FontSize', 14);
ylabel('Y坐标 (m)', 'FontSize', 14);

% 绘制节点
scatter(xCoords, yCoords, 80, 'k', 'filled', 'MarkerEdgeColor', 'w', 'LineWidth', 1.5);

% 添加节点标签
for i = 1:numNodes
    text(xCoords(i)+3, yCoords(i)+3, sprintf('%d', nodeIDs(i)), ...
         'FontSize', 9, 'FontWeight', 'bold', 'Color', 'w', ...
         'BackgroundColor', [0,0,0,0.5], ...
         'HorizontalAlignment', 'center');
end

% 标记放射源
scatter(sources(:,1), sources(:,2), 150, 'r', 'pentagram', 'filled', 'LineWidth', 1.5);

% 添加等高线
contour(xGrid, yGrid, radGrid, 10, 'LineColor', 'k', 'LineWidth', 0.5);

% 添加图例
legend('节点位置', '放射源位置', 'Location', 'northeast', 'FontSize', 12);

% 设置坐标轴
axis equal;
xlim([0 200]);
ylim([0 200]);

hold off;

%% 输出结果
fprintf('\n===== 案例2复现结果 (DBO算法) =====\n');
fprintf('最优路径: %s\n', strjoin(cellstr(num2str(gBestPath')), '->'));
fprintf('最优路径总暴露代价: %.2f a.u.\n', gBestFitness);
fprintf('算法总迭代次数: %d\n', maxIter);

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