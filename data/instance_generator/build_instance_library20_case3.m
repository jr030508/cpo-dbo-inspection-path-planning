function build_instance_library20_case3()
% build_instance_library20_case3.m
% ============================================================
%  Case3 (51 nodes) 20个扰动实例库母本（统一场景，供所有算法对照）
%  输出（默认输出到 ./instance_lib20/ ）：
%   1) instance_library20.mat  (含 doseMatrices: 51×51×20)
%   2) doseMatrix_inst01.csv ... doseMatrix_inst20.csv
%   3) instance_summary_20.csv  (每个实例的 noise/strength/sources 摘要)
%
%  扰动规则（每个实例都有）：
%   (1) 放射源强度 ±10% 扰动，并整体缩放保持 sum(strength) 不变（均值不变）
%   (2) 放射源位置 ±posJitter (默认 2m) 均匀扰动，并裁剪到[0,200]
%   (3) doseMatrix 加相对误差噪声：1%~3%（对称、非负、对角线0）
%   (4) 热点增强：top 5% 边剂量 ×3（与之前保持一致）
%
%  复现性：
%   rng(6,'twister') 固定，生成的20个实例永远一致
%
%  重要修复（无坑版）：
%   ✅ 1) 文件函数化：避免脚本末尾 end 报错
%   ✅ 2) instSources 摘要 x/y 顺序修复：x1..x7 与 y1..y7 不再错位
% ============================================================

clc; close all;
rng(6,'twister');

%% ============ 0) 实例设置 ============
numInstances    = 20;
posJitter       = 2;      % ±2m（想更强可改 5）
strengthJitter  = 0.10;   % ±10%
noiseMin        = 0.01;   % 1%
noiseMax        = 0.03;   % 3%

% 输出目录（建议集中管理，避免文件散落）
outDir = fullfile(pwd, 'instance_lib20');
if ~exist(outDir), mkdir(outDir); end

outMat = fullfile(outDir, 'instance_library20.mat');

fprintf('\n================= 生成 20 个扰动实例库（rng=6）=================\n');
fprintf('Output folder: %s\n', outDir);

%% ============ 1) 51个节点 ============
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

%% ============ 2) 基准放射源 ============
sources0 = [
    35, 35;
    30, 110;
    88, 70;
    88, 142;
    138, 107;
    160, 50;
    168, 180;
];
strength0 = [8000, 12000, 10000, 15000, 9000, 13000, 11000];

assert(size(sources0,1)==numel(strength0), 'sources0 行数与 strength0 长度不一致。');
assert(size(sources0,2)==2, 'sources0 必须是 N×2。');

%% ============ 3) 预分配（摘要用 x/y 分开存，彻底避免错位） ============
instNoiseLevel = zeros(numInstances,1);
instStrengths  = zeros(numInstances, numel(strength0));
instX          = zeros(numInstances, size(sources0,1));  % 20×7
instY          = zeros(numInstances, size(sources0,1));  % 20×7

doseMatrices   = zeros(numNodes, numNodes, numInstances);

%% ============ 4) 生成 20 个实例 ============
for k = 1:numInstances
    [sources_k, strength_k, noise_k] = perturbOneInstance( ...
        sources0, strength0, posJitter, strengthJitter, noiseMin, noiseMax);

    instNoiseLevel(k)   = noise_k;
    instStrengths(k,:)  = strength_k(:)';

    % ✅ 修复点：x/y 分列存储，instance_summary.csv 不会再错
    instX(k,:) = sources_k(:,1)';
    instY(k,:) = sources_k(:,2)';

    dose_k = buildDoseMatrix_case3(xCoords, yCoords, sources_k, strength_k, noise_k);
    doseMatrices(:,:,k) = dose_k;

    csvName = fullfile(outDir, sprintf('doseMatrix_inst%02d.csv', k));
    writematrix(dose_k, csvName);

    fprintf('Inst %02d/%02d | noise=%.3f | sumStrength=%.1f | saved=%s\n', ...
        k, numInstances, noise_k, sum(strength_k), getFileName(csvName));
end

%% ============ 5) 导出摘要表（列名严格对应） ============
% Variable names：Instance, NoiseLevel, S1..S7, x1..x7, y1..y7
varNames = [{'Instance','NoiseLevel'}, ...
            arrayfun(@(i)sprintf('S%d',i),1:7,'UniformOutput',false), ...
            arrayfun(@(i)sprintf('x%d',i),1:7,'UniformOutput',false), ...
            arrayfun(@(i)sprintf('y%d',i),1:7,'UniformOutput',false)];

T = table((1:numInstances)', instNoiseLevel, ...
          instStrengths(:,1),instStrengths(:,2),instStrengths(:,3),instStrengths(:,4),instStrengths(:,5),instStrengths(:,6),instStrengths(:,7), ...
          instX(:,1),instX(:,2),instX(:,3),instX(:,4),instX(:,5),instX(:,6),instX(:,7), ...
          instY(:,1),instY(:,2),instY(:,3),instY(:,4),instY(:,5),instY(:,6),instY(:,7), ...
          'VariableNames', varNames);

summaryCsv = fullfile(outDir, 'instance_summary_20.csv');
writetable(T, summaryCsv);

%% ============ 6) 保存 MAT（核心：包含 3D doseMatrices） ============
save(outMat, ...
    'numInstances','posJitter','strengthJitter','noiseMin','noiseMax', ...
    'nodeData','xCoords','yCoords','sources0','strength0', ...
    'instNoiseLevel','instStrengths','instX','instY','doseMatrices');

fprintf('\n✅ 母本已生成并保存：%s\n', outMat);
fprintf('✅ doseMatrices 尺寸：%s（应为 51x51x20）\n', mat2str(size(doseMatrices)));
fprintf('✅ 同时导出：doseMatrix_inst01..20.csv + instance_summary_20.csv\n');
fprintf('✅ 摘要表列顺序：Instance, NoiseLevel, S1..S7, x1..x7, y1..y7（已对齐）\n');

end % ===== end main function =====


%% ====================== Local functions ======================

function [sources, strength, noiseLevel] = perturbOneInstance(sources0, strength0, posJitter, strengthJitter, noiseMin, noiseMax)
% (1) 强度扰动 ±strengthJitter，然后缩放保持 sum(strength) 不变
m = numel(strength0);
f = 1 + (rand(1,m)*2 - 1) * strengthJitter;
strength_raw = strength0(:)' .* f;

% 防止极端情况（虽然这里不会出现负数，但做个保护）
strength_raw(strength_raw < eps) = eps;

strength = strength_raw * (sum(strength0) / sum(strength_raw));

% (2) 位置扰动 ±posJitter，并裁剪到[0,200]
dxy = (rand(size(sources0))*2 - 1) * posJitter;
sources = sources0 + dxy;
sources(:,1) = min(200, max(0, sources(:,1)));
sources(:,2) = min(200, max(0, sources(:,2)));

% (3) 噪声水平 noiseMin~noiseMax
noiseLevel = noiseMin + (noiseMax-noiseMin)*rand();
end


function doseMatrix = buildDoseMatrix_case3(xCoords, yCoords, sources, sourceStrength, noiseLevel)
% 构建 doseMatrix + 对称相对噪声 + top5%×3 热点增强

numNodes = numel(xCoords);
doseMatrix = zeros(numNodes, numNodes);

% 1) 节点辐射水平
nodeRadLevel = zeros(numNodes, 1);
for i = 1:numNodes
    totalRad = 0;
    for s = 1:size(sources, 1)
        dist = hypot(xCoords(i)-sources(s,1), yCoords(i)-sources(s,2));
        totalRad = totalRad + sourceStrength(s) / (dist^2 + 10);
    end
    nodeRadLevel(i) = totalRad;
end

% 2) 边剂量：端点均值 + 中点项
for i = 1:numNodes
    for j = i+1:numNodes
        dist = hypot(xCoords(i)-xCoords(j), yCoords(i)-yCoords(j));
        avgRad = (nodeRadLevel(i) + nodeRadLevel(j)) / 2;

        midX = (xCoords(i) + xCoords(j)) / 2;
        midY = (yCoords(i) + yCoords(j)) / 2;
        midRad = 0;
        for s = 1:size(sources, 1)
            distToSource = hypot(midX-sources(s,1), midY-sources(s,2));
            midRad = midRad + (sourceStrength(s)/2) / (distToSource^2 + 10);
        end

        dose = dist * (avgRad + midRad) / 100;
        doseMatrix(i,j) = dose;
        doseMatrix(j,i) = dose;
    end
end
doseMatrix(1:numNodes+1:end) = 0;

% 3) 相对噪声（对称、非负、对角线0）
E = (rand(numNodes) * 2 - 1) * noiseLevel;  % [-p,+p]
E = triu(E,1);
E = E + E.';                                % 对称
doseMatrix = doseMatrix .* (1 + E);
doseMatrix(doseMatrix < 0) = 0;
doseMatrix(1:numNodes+1:end) = 0;

% 4) 热点增强 top 5% ×3（只对非零边）
allD = doseMatrix(:);
allD(allD == 0) = [];
if ~isempty(allD)
    thr = prctile(allD, 95);
    for i = 1:numNodes
        for j = i+1:numNodes
            if doseMatrix(i,j) > thr
                doseMatrix(i,j) = doseMatrix(i,j) * 3;
                doseMatrix(j,i) = doseMatrix(j,i) * 3;
            end
        end
    end
end
doseMatrix(1:numNodes+1:end) = 0;

end


function name = getFileName(fullpath)
% 只用于打印：输出文件名而非整段路径
[~, name, ext] = fileparts(fullpath);
name = [name ext];
end