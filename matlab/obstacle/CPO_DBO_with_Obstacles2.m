%% ============================================================
%  Case3 (51 nodes) 核场景巡检路径规划（含障碍） - 混合边模型 + 新CPO-DBO(你给的版本)
%
%  混合边模型（你的思路）：
%    - 若两点之间无遮挡（LOS，相对扩张障碍矩形）：
%        -> 直线边（效率）
%        -> 边剂量 = 沿直线对辐射场采样积分（trapz）
%    - 若被障碍遮挡：
%        -> 栅格 A* 搜索可行路径（不穿障）
%        -> 目标：剂量最小 + epsLen*距离（极小 tie-break，避免奇怪绕路）
%
%  优化算法：CPO-DBO hybrid（Permutation-safe，采用你提供的“新CPO-DBO”实现）
%    - OX交叉 + swap/insert/inversion(+scramble用于escape)
%    - 动态概率 + 分组偏置 + 精英保留 + 轻量2-opt
%
%  输出文件：
%    - cost_matrix_mixed_case3.csv      (Dose + epsLen*Len)
%    - dose_matrix_mixed_case3.csv
%    - pathlen_matrix_mixed_case3.csv
%    - edgeType_matrix_mixed_case3.csv (0=LOS直线, 1=A*绕障)
%    - obstacles_case3.mat
%
%  出图（中文 + 图例修复 + 场景图 + 剂量分解条形图）：
%    1) 收敛曲线
%    2) 场景图（节点+点源+障碍，不含路径）
%    3) 最优路径图（LOS实线 / A*虚线，图例文字完整）
%    4) 剂量分解条形图（逐边剂量贡献 + 标红A*边）
%% ============================================================

clc; clear; close all;
rng(6);

%% ===================== 0) USER PARAMETERS ===================
% ---- Primary objective: DOSE ----
% epsLen is a tiny tie-break for length (keeps A* from making silly long loops)
epsLen = 1e-5;

% ---- Grid resolution (m per cell) ----
gridRes = 2.0;

% ---- Straight-line dose integration step (meters) ----
dsLine = 1.0;

% ---- CPO-DBO parameters (use your "new CPO-DBO" style) ----
M = 400;
maxIter = 5000;
%Tscale = 30000;                 % 固定节奏：以30000为标尺


eliteKeep = max(2, round(0.057*M));   % 精英保留
twoOptIters_defense = 1;
twoOptIters_attack  = 1;
twoOptIters_escape  = 0;

% DBO group split
popBall  = round(M * 0.4);
popBreed = round(M * 0.2);
popBaby  = M - popBall - popBreed;

% dynamic behavior probabilities
defenseProb0 = 0.70;
attackProb0  = 0.25;
escapeProb0  = 0.05;

% ---- Obstacles ----
% format: [xmin, ymin, width, height]
obstacles = [
    55,  12,  14,  10;   % 障碍1
    96,  22,  18,  14;   % 障碍2
   % 92,  40,  18,  14;   % 障碍2
    5,   115, 14,  14;   % 障碍3
    110, 150, 16,  18;   % 障碍4
    150, 120, 18,  14;   % 障碍5
   % 150, 120, 18,  14;   % 障碍5
];
obsMargin = 2.0; % safety buffer (expanded obstacles used for LOS + occupancy)

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

%% ===================== 2) Obstacle overlap check & save ======
assertObstaclesNoOverlap(obstacles, obsMargin, xCoords, yCoords, sources);
save('obstacles_case3.mat', 'obstacles', 'obsMargin');

%% ===================== 3) Build grid, occupancy, radiation ===
xMin = 0; xMax = 200;
yMin = 0; yMax = 200;

xVec = xMin:gridRes:xMax;
yVec = yMin:gridRes:yMax;
Ny = numel(yVec);
Nx = numel(xVec);

occ = false(Ny, Nx); % row=y, col=x
occ = rasterizeObstacles(occ, obstacles, obsMargin, xVec, yVec);

[XX, YY] = meshgrid(xVec, yVec);
radGrid = zeros(Ny, Nx);
for s = 1:size(sources,1)
    dx = XX - sources(s,1);
    dy = YY - sources(s,2);
    radGrid = radGrid + sourceStrength(s) ./ (dx.^2 + dy.^2 + 10);
end
minRad = min(radGrid(:));

% map nodes to grid indices
nodeRC = worldToGrid(xCoords, yCoords, xVec, yVec);
for i = 1:numNodes
    if occ(nodeRC(i,1), nodeRC(i,2))
        error('节点 %d 映射到了障碍占用栅格内，请调整 obstacles/gridRes。', i);
    end
end

%% ===================== 4) Build MIXED pairwise matrices ======
doseMat  = inf(numNodes, numNodes);
lenMat   = inf(numNodes, numNodes);
costMat  = inf(numNodes, numNodes);
edgeType = zeros(numNodes, numNodes); % 0=LOS straight, 1=A*

fprintf('\n[STEP] 构建混合边矩阵：无遮挡->直线，遮挡->A*(剂量最小)...\n');
totalPairs = numNodes*(numNodes-1)/2;
cnt = 0;

for i = 1:numNodes
    for j = i+1:numNodes
        cnt = cnt + 1;

        p1 = [xCoords(i), yCoords(i)];
        p2 = [xCoords(j), yCoords(j)];

        isLOS = hasLOS_world(p1, p2, obstacles, obsMargin);

        if isLOS
            L = hypot(p2(1)-p1(1), p2(2)-p1(2));
            D = lineDoseIntegral(p1, p2, sources, sourceStrength, dsLine);
            edgeType(i,j) = 0; edgeType(j,i) = 0;
        else
            [~, L, D] = astar_grid_dose(nodeRC(i,:), nodeRC(j,:), occ, radGrid, gridRes, epsLen, minRad);
            edgeType(i,j) = 1; edgeType(j,i) = 1;
        end

        lenMat(i,j)  = L; lenMat(j,i)  = L;
        doseMat(i,j) = D; doseMat(j,i) = D;

        costMat(i,j) = D + epsLen * L;
        costMat(j,i) = costMat(i,j);

        if mod(cnt,100)==0
            fprintf('  Mixed edges: %d/%d pairs done\n', cnt, totalPairs);
        end
    end
end

for i=1:numNodes
    lenMat(i,i)=0; doseMat(i,i)=0; costMat(i,i)=0; edgeType(i,i)=0;
end

writematrix(lenMat,  'pathlen_matrix_mixed_case3.csv');
writematrix(doseMat, 'dose_matrix_mixed_case3.csv');
writematrix(costMat, 'cost_matrix_mixed_case3.csv');
writematrix(edgeType,'edgeType_matrix_mixed_case3.csv');

fprintf('\nSaved:\n- pathlen_matrix_mixed_case3.csv\n- dose_matrix_mixed_case3.csv\n- cost_matrix_mixed_case3.csv\n- edgeType_matrix_mixed_case3.csv\n');
fprintf('Settings: epsLen=%.2e, gridRes=%.2f, dsLine=%.2f, obsMargin=%.2f\n', epsLen, gridRes, dsLine, obsMargin);

%% ===================== 5) 新 CPO-DBO on mixed cost matrix =====
population = cell(M,1);
fitness    = zeros(M,1);

for i = 1:M
    perm = randperm(numNodes-1)+1;
    population{i} = [1, perm, 1];
    fitness(i) = tourCost(population{i}, costMat);
end

pBestPath = population;
pBestFitness = fitness;
[gBestFitness, gBestIdx] = min(fitness);
gBestPath = population{gBestIdx};

convergence = zeros(maxIter,1);

fprintf('\n[STEP] Start 新CPO-DBO on MIXED cost matrix (Dose + epsLen*Len)...\n');

for iter = 1:maxIter
    t = iter / maxIter;
    %t = min(1, iter / Tscale);      % 5000代时 t 只到 0.166

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

        % 选择行为（按你给的逻辑）
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

%% ===================== 6) Output summary =====================
totalDose = tourCost(gBestPath, doseMat);
totalLen  = tourCost(gBestPath, lenMat);
nAstarEdges = countAstarEdgesInTour(gBestPath, edgeType);

fprintf('\n===== Case3 混合边（LOS直线 + A*绕障）+ 新CPO-DBO =====\n');
fprintf('epsLen=%.2e (tie-break), gridRes=%.2f, dsLine=%.2f\n', epsLen, gridRes, dsLine);
fprintf('Best tour COST (Dose + epsLen*Len): %.6f\n', gBestFitness);
fprintf('Total Dose (main): %.6f\n', totalDose);
fprintf('Total Length (ref): %.6f\n', totalLen);
fprintf('A* edges in best tour: %d / %d\n', nAstarEdges, numNodes);
fprintf('Best tour: %s\n', strjoin(cellstr(num2str(gBestPath')), '->'));

%% ===================== 7) Visualization (中文) =================
% 7.1 收敛曲线
figure('Name','收敛曲线（混合边+新CPO-DBO）','Position',[100,100,900,420]);
plot(convergence,'LineWidth',2); grid on;
xlabel('迭代次数');
ylabel('最优巡检代价（剂量 + \epsilon \cdot 距离）');
title('新CPO-DBO 在混合边权矩阵上的收敛曲线');

% 7.2 场景图（不含路径）
figure('Name','场景示意图（不含路径）','Position',[100,100,1000,820]);
hold on; axis equal; grid on;
xlim([xMin xMax]); ylim([yMin yMax]);
xlabel('X坐标'); ylabel('Y坐标');
title('核场景地图：节点、点源与障碍物（不含路径）');

contourf(XX, YY, radGrid, 30, 'LineStyle','none');
colormap(jet); cb = colorbar; cb.Label.String = '辐射水平（相对值）';

drawObstacles_CN(obstacles);

scatter(xCoords, yCoords, 70, 'k', 'filled', 'MarkerEdgeColor','w');
for i = 1:numNodes
    text(xCoords(i)+2, yCoords(i)+1, sprintf('%d',i), ...
        'Color','w','FontWeight','bold', 'BackgroundColor',[0 0 0 0.35], ...
        'HorizontalAlignment','center');
end

scatter(sources(:,1), sources(:,2), 160, 'r', 'pentagram', 'filled');
for s=1:size(sources,1)
    text(sources(s,1), sources(s,2)-5, sprintf('源%d',s), ...
        'Color','r','FontWeight','bold', 'HorizontalAlignment','center');
end

hLegObs  = plot(nan,nan,'s','MarkerFaceColor',[0.3 0.3 0.3],'MarkerEdgeColor','k','LineStyle','none');
hLegNode = plot(nan,nan,'o','MarkerFaceColor','k','MarkerEdgeColor','w','LineStyle','none');
hLegSrc  = plot(nan,nan,'p','MarkerFaceColor','r','MarkerEdgeColor','r','LineStyle','none');
legend([hLegObs, hLegNode, hLegSrc], {'障碍物','巡检节点','点源'}, 'Location','northeastoutside');
hold off;

% 7.3 最优路径图（图例修复）
figure('Name','最优巡检路径（混合边+新CPO-DBO）','Position',[100,100,1150,880]);
hold on; axis equal; grid on;
xlim([xMin xMax]); ylim([yMin yMax]);
xlabel('X坐标'); ylabel('Y坐标');
title(sprintf('最优巡检回路（总剂量=%.3f，总距离=%.1f，A*边=%d，\\epsilon=%.1e）', ...
    totalDose, totalLen, nAstarEdges, epsLen));

contourf(XX, YY, radGrid, 30, 'LineStyle','none');
colormap(jet); cb = colorbar; cb.Label.String = '辐射水平（相对值）';

drawObstacles_CN(obstacles);

scatter(xCoords, yCoords, 70, 'k', 'filled', 'MarkerEdgeColor','w');
for i = 1:numNodes
    text(xCoords(i)+2, yCoords(i)+1, sprintf('%d',i), ...
        'Color','w','FontWeight','bold', 'BackgroundColor',[0 0 0 0.35], ...
        'HorizontalAlignment','center');
end

scatter(sources(:,1), sources(:,2), 160, 'r', 'pentagram', 'filled');
for s=1:size(sources,1)
    text(sources(s,1), sources(s,2)-5, sprintf('源%d',s), ...
        'Color','r','FontWeight','bold', 'HorizontalAlignment','center');
end

cols = lines(length(gBestPath)-1);
for k = 1:(length(gBestPath)-1)
    a = gBestPath(k);
    b = gBestPath(k+1);

    p1 = [xCoords(a), yCoords(a)];
    p2 = [xCoords(b), yCoords(b)];

    if edgeType(a,b) == 0
        plot([p1(1) p2(1)], [p1(2) p2(2)], ...
            'LineWidth', 2.2, 'LineStyle','-', ...
            'Color', cols(mod(k-1,size(cols,1))+1,:));
    else
        [pathRC,~,~] = astar_grid_dose(nodeRC(a,:), nodeRC(b,:), occ, radGrid, gridRes, epsLen, minRad);
        pathRC2 = simplifyCollinear(pathRC);
        pathRC3 = shortcutLineOfSight(pathRC2, occ);
        [px, py] = gridPathToWorld(pathRC3, xVec, yVec);
        plot(px, py, ...
            'LineWidth', 2.4, 'LineStyle','--', ...
            'Color', cols(mod(k-1,size(cols,1))+1,:));
    end
end

hLegObs  = plot(nan,nan,'s','MarkerFaceColor',[0.3 0.3 0.3],'MarkerEdgeColor','k','LineStyle','none');
hLegNode = plot(nan,nan,'o','MarkerFaceColor','k','MarkerEdgeColor','w','LineStyle','none');
hLegSrc  = plot(nan,nan,'p','MarkerFaceColor','r','MarkerEdgeColor','r','LineStyle','none');
hLegLOS  = plot(nan,nan,'-','LineWidth',2.2,'Color','k');
hLegASt  = plot(nan,nan,'--','LineWidth',2.4,'Color','k');

legend([hLegObs,hLegNode,hLegSrc,hLegLOS,hLegASt], ...
    {'障碍物','巡检节点','点源','无遮挡直线边（LOS）','绕障边（A*）'}, ...
    'Location','northeastoutside');
hold off;

% 7.4 剂量分解条形图
[edgeDoseList, edgeTypeList] = getDoseBreakdown(gBestPath, doseMat, edgeType);
figure('Name','剂量分解条形图','Position',[120,120,1100,420]);
bar(edgeDoseList);
grid on;
xlabel('回路边序号');
ylabel('该边剂量贡献');
title('最优巡检回路的剂量分解（逐边贡献）');
hold on;
idxA = find(edgeTypeList==1);
if ~isempty(idxA)
    scatter(idxA, edgeDoseList(idxA), 50, 'r', 'filled', 'MarkerEdgeColor','k');
end
hold off;
legend({'边剂量','A*绕障边'}, 'Location','northeast');

%% ===================== 8) Save MAT result =====================
resultFile = 'result_case3_mixed_obstacle.mat';

save(resultFile, ...
    'gBestPath', 'gBestFitness', 'convergence', ...
    'doseMat', 'lenMat', 'costMat', 'edgeType', ...
    'xCoords', 'yCoords', 'nodeData', 'numNodes', ...
    'sources', 'sourceStrength', ...
    'obstacles', 'obsMargin', ...
    'epsLen', 'gridRes', 'dsLine', ...
    'xMin', 'xMax', 'yMin', 'yMax', ...
    'xVec', 'yVec', 'XX', 'YY', ...
    'occ', 'nodeRC', 'radGrid', 'minRad', ...
    'M', 'maxIter', 'eliteKeep', ...
    'twoOptIters_defense', 'twoOptIters_attack', 'twoOptIters_escape', ...
    'popBall', 'popBreed', 'popBaby', ...
    'defenseProb0', 'attackProb0', 'escapeProb0', ...
    'totalDose', 'totalLen', 'nAstarEdges', ...
    'edgeDoseList', 'edgeTypeList');

fprintf('\n已保存 MAT 结果文件: %s\n', resultFile);
%                      Functions
%% ============================================================

function assertObstaclesNoOverlap(obstacles, margin, xNodes, yNodes, sources)
    ptsX = [xNodes(:); sources(:,1)];
    ptsY = [yNodes(:); sources(:,2)];
    labels = [repmat("节点", numel(xNodes), 1); repmat("点源", size(sources,1), 1)];
    for k = 1:size(obstacles,1)
        rect = obstacles(k,:);
        rectExp = [rect(1)-margin, rect(2)-margin, rect(3)+2*margin, rect(4)+2*margin];
        inside = pointInRect(ptsX, ptsY, rectExp);
        if any(inside)
            idx = find(inside, 1, 'first');
            error('障碍物%d 与%s重合（x=%.1f,y=%.1f）。请调整 obstacles 或 obsMargin。', ...
                k, labels(idx), ptsX(idx), ptsY(idx));
        end
    end
    fprintf('\n障碍检查通过：障碍物未与节点/点源重合（含margin=%.1f）\n', margin);
end

function inside = pointInRect(px, py, rect)
    xmin = rect(1); ymin = rect(2); w = rect(3); h = rect(4);
    inside = (px >= xmin) & (px <= xmin+w) & (py >= ymin) & (py <= ymin+h);
end

function occ = rasterizeObstacles(occ, obstacles, margin, xVec, yVec)
    [XX, YY] = meshgrid(xVec, yVec);
    for k = 1:size(obstacles,1)
        r = obstacles(k,:);
        xmin = r(1)-margin; ymin = r(2)-margin;
        xmax = r(1)+r(3)+margin; ymax = r(2)+r(4)+margin;
        mask = (XX>=xmin) & (XX<=xmax) & (YY>=ymin) & (YY<=ymax);
        occ(mask) = true;
    end
end

function drawObstacles_CN(obstacles)
    for k = 1:size(obstacles,1)
        r = obstacles(k,:);
        rectangle('Position', r, 'FaceColor', [0.3 0.3 0.3 0.65], ...
            'EdgeColor','k', 'LineWidth',1.5);
        cx = r(1)+r(3)/2; cy = r(2)+r(4)/2;
        text(cx, cy, sprintf('障碍%d',k), 'Color','w','FontWeight','bold', ...
            'HorizontalAlignment','center','VerticalAlignment','middle');
    end
end

function rc = worldToGrid(x, y, xVec, yVec)
    rc = zeros(numel(x),2);
    for i = 1:numel(x)
        [~, c] = min(abs(xVec - x(i)));
        [~, r] = min(abs(yVec - y(i)));
        rc(i,:) = [r, c];
    end
end

function [px, py] = gridPathToWorld(pathRC, xVec, yVec)
    py = yVec(pathRC(:,1));
    px = xVec(pathRC(:,2));
end

function c = tourCost(path, mat)
    c = 0;
    for i = 1:(length(path)-1)
        c = c + mat(path(i), path(i+1));
    end
end

function n = countAstarEdgesInTour(tour, edgeType)
    n = 0;
    for k=1:(numel(tour)-1)
        if edgeType(tour(k), tour(k+1)) == 1
            n = n + 1;
        end
    end
end

function [doseList, typeList] = getDoseBreakdown(tour, doseMat, edgeType)
    m = numel(tour)-1;
    doseList = zeros(m,1);
    typeList = zeros(m,1);
    for k=1:m
        a=tour(k); b=tour(k+1);
        doseList(k)=doseMat(a,b);
        typeList(k)=edgeType(a,b);
    end
end

%% ---------- LOS check in continuous geometry ----------
function ok = hasLOS_world(p1, p2, obstacles, margin)
    ok = true;
    for k=1:size(obstacles,1)
        r = obstacles(k,:);
        rect = [r(1)-margin, r(2)-margin, r(3)+2*margin, r(4)+2*margin];
        if segmentIntersectsRect(p1, p2, rect)
            ok = false; return;
        end
    end
end

function hit = segmentIntersectsRect(p1, p2, rect)
    xmin=rect(1); ymin=rect(2); xmax=xmin+rect(3); ymax=ymin+rect(4);

    segxmin=min(p1(1),p2(1)); segxmax=max(p1(1),p2(1));
    segymin=min(p1(2),p2(2)); segymax=max(p1(2),p2(2));
    if segxmax < xmin || segxmin > xmax || segymax < ymin || segymin > ymax
        hit = false; return;
    end

    if inRect(p1,xmin,ymin,xmax,ymax) || inRect(p2,xmin,ymin,xmax,ymax)
        hit = true; return;
    end

    A=[xmin ymin]; B=[xmax ymin];
    C=[xmax ymax]; D=[xmin ymax];
    hit = segmentsIntersect(p1,p2,A,B) || segmentsIntersect(p1,p2,B,C) || ...
          segmentsIntersect(p1,p2,C,D) || segmentsIntersect(p1,p2,D,A);
end

function tf = inRect(p, xmin,ymin,xmax,ymax)
    tf = (p(1)>=xmin && p(1)<=xmax && p(2)>=ymin && p(2)<=ymax);
end

function tf = segmentsIntersect(p1,p2,p3,p4)
    tf = (ccw(p1,p3,p4) ~= ccw(p2,p3,p4)) && (ccw(p1,p2,p3) ~= ccw(p1,p2,p4));
end

function v = ccw(a,b,c)
    v = (c(2)-a(2))*(b(1)-a(1)) > (b(2)-a(2))*(c(1)-a(1));
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
    for s=1:size(sources,1)
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

%% ---------- A* on grid: MINIMUM DOSE (with epsLen tie-break) ----------
function [pathRC, pathLen, pathDose] = astar_grid_dose(startRC, goalRC, occ, radGrid, gridRes, epsLen, minRad)

    Ny = size(occ,1); Nx = size(occ,2); %#ok<NASGU>

    moves = [
        -1,  0, 1;
         1,  0, 1;
         0, -1, 1;
         0,  1, 1;
        -1, -1, sqrt(2);
        -1,  1, sqrt(2);
         1, -1, sqrt(2);
         1,  1, sqrt(2);
    ];

    sr = startRC(1); sc = startRC(2);
    gr = goalRC(1);  gc = goalRC(2);

    gScore = inf(size(occ));
    fScore = inf(size(occ));
    cameR = zeros(size(occ), 'int16');
    cameC = zeros(size(occ), 'int16');
    inOpen = false(size(occ));
    closed = false(size(occ));

    gScore(sr,sc) = 0;

    hPerMeter = (minRad/100) + epsLen;
    fScore(sr,sc) = hPerMeter * (gridRes * hypot(gr-sr, gc-sc));
    inOpen(sr,sc) = true;

    while any(inOpen(:))
        [cr, cc] = argminOpen(fScore, inOpen);
        if cr==gr && cc==gc
            break;
        end

        inOpen(cr,cc) = false;
        closed(cr,cc) = true;

        for k = 1:size(moves,1)
            nr = cr + moves(k,1);
            nc = cc + moves(k,2);
            if nr<1 || nr>size(occ,1) || nc<1 || nc>size(occ,2)
                continue;
            end
            if occ(nr,nc) || closed(nr,nc)
                continue;
            end

            stepLen  = moves(k,3) * gridRes;
            stepDose = stepLen * radGrid(nr,nc) / 100;
            stepCost = stepDose + epsLen * stepLen;

            tentativeG = gScore(cr,cc) + stepCost;

            if tentativeG < gScore(nr,nc)
                cameR(nr,nc) = cr;
                cameC(nr,nc) = cc;
                gScore(nr,nc) = tentativeG;

                h = hPerMeter * (gridRes * hypot(gr-nr, gc-nc));
                fScore(nr,nc) = tentativeG + h;
                inOpen(nr,nc) = true;
            end
        end
    end

    if ~isfinite(gScore(gr,gc))
        error('A* failed: (%d,%d)->(%d,%d) unreachable.', sr,sc,gr,gc);
    end

    pathRC = reconstructPath(cameR, cameC, startRC, goalRC);

    pathLen = 0; pathDose = 0;
    for i = 1:(size(pathRC,1)-1)
        dr = abs(pathRC(i+1,1)-pathRC(i,1));
        dc = abs(pathRC(i+1,2)-pathRC(i,2));
        if dr==1 && dc==1
            stepLen = sqrt(2)*gridRes;
        else
            stepLen = 1*gridRes;
        end
        pathLen = pathLen + stepLen;
        r = pathRC(i+1,1); c = pathRC(i+1,2);
        pathDose = pathDose + stepLen * radGrid(r,c) / 100;
    end
end

function [r,c] = argminOpen(fScore, inOpen)
    tmp = fScore; tmp(~inOpen) = inf;
    [~, idx] = min(tmp(:));
    [r,c] = ind2sub(size(tmp), idx);
end

function pathRC = reconstructPath(cameR, cameC, startRC, goalRC)
    cur = goalRC;
    pathRC = cur;
    while ~(cur(1)==startRC(1) && cur(2)==startRC(2))
        pr = cameR(cur(1),cur(2));
        pc = cameC(cur(1),cur(2));
        if pr==0 && pc==0
            error('Path reconstruction failed.');
        end
        cur = [double(pr), double(pc)];
        pathRC = [cur; pathRC]; %#ok<AGROW>
    end
end

%% ---------- Plot smoothing (only for visualization) ----------
function path2 = simplifyCollinear(pathRC)
    if size(pathRC,1) <= 2
        path2 = pathRC; return;
    end
    keep = true(size(pathRC,1),1);
    for i = 2:size(pathRC,1)-1
        v1 = pathRC(i,:) - pathRC(i-1,:);
        v2 = pathRC(i+1,:) - pathRC(i,:);
        if v1(1)==v2(1) && v1(2)==v2(2)
            keep(i) = false;
        end
    end
    path2 = pathRC(keep,:);
end

function path3 = shortcutLineOfSight(pathRC, occ)
    if size(pathRC,1) <= 2
        path3 = pathRC; return;
    end
    path3 = pathRC(1,:);
    i = 1;
    while i < size(pathRC,1)
        j = size(pathRC,1);
        while j > i+1
            if hasLineOfSight_grid(pathRC(i,:), pathRC(j,:), occ)
                break;
            end
            j = j - 1;
        end
        path3 = [path3; pathRC(j,:)]; %#ok<AGROW>
        i = j;
    end
end

function ok = hasLineOfSight_grid(aRC, bRC, occ)
    [rr, cc] = bresenham(aRC(1), aRC(2), bRC(1), bRC(2));
    ok = true;
    for k = 1:numel(rr)
        if occ(rr(k), cc(k))
            ok = false; return;
        end
    end
end

function [rr, cc] = bresenham(r0, c0, r1, c1)
    dr = abs(r1 - r0);
    dc = abs(c1 - c0);
    sr = sign(r1 - r0); if sr==0, sr=1; end
    sc = sign(c1 - c0); if sc==0, sc=1; end

    if dc > dr
        [cc, rr] = bresenhamCore(c0, r0, c1, r1, sc, sr, dc, dr);
    else
        [rr, cc] = bresenhamCore(r0, c0, r1, c1, sr, sc, dr, dc);
    end
end

function [aa, bb] = bresenhamCore(a0, b0, a1, b1, sa, sb, da, db)
    aa = zeros(da+1,1);
    bb = zeros(da+1,1);
    err = 2*db - da;
    a = a0; b = b0;
    for i = 1:(da+1)
        aa(i) = a; bb(i) = b;
        if err >= 0
            b = b + sb;
            err = err - 2*da;
        end
        a = a + sa;
        err = err + 2*db;
    end
end