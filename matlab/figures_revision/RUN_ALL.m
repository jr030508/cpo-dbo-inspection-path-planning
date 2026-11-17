%% ============================================================
%  一键重绘论文 7 张需要修改标注的图片
%  使用方法：在 MATLAB 中打开本文件，点击上方"运行"按钮即可
%  全部图片输出到 output 文件夹（PNG，300dpi，Times New Roman 字体）
%% ============================================================
clc; clear; close all;
scriptDir = fileparts(mfilename('fullpath'));
cd(scriptDir);
fprintf('=============================================================\n');
fprintf(' 论文图片重绘（相对暴露代价 a.u. 版，Times New Roman）\n');
fprintf(' 工作目录: %s\n', scriptDir);
fprintf('=============================================================\n\n');

jobs = {'fig4_heatmap_cost.m',        '图4  30节点暴露代价矩阵热力图';
        'fig6_convergence_cost.m',    '图6  四算法收敛曲线';
        'fig8_boxplots_cost.m',       '图8  20实例四指标箱线图';
        'fig10_ecdf_cost.m',          '图10 ECDF曲线';
        'fig11_risk_bars_cost.m',     '图11 P90/P95/Max风险柱状图';
        'fig14_seed_curves_cost.m',   '图14 五种子收敛曲线';
        'fig21_free_vs_obstacle_cost.m', '图21 无障碍vs有障碍路径对比'};

nOK = 0;
for k = 1:size(jobs,1)
    fprintf('[%d/7] 正在绘制 %s ...\n', k, jobs{k,2});
    try
        run(fullfile(scriptDir, jobs{k,1}));
        nOK = nOK + 1;
    catch ME
        fprintf('   !! 失败: %s\n   !! 原因: %s\n', jobs{k,1}, ME.message);
    end
end

fprintf('\n=============================================================\n');
fprintf(' 全部完成: %d / 7 张图输出到 output 文件夹\n', nOK);
fprintf(' 请将 PNG 按对照表插回 Word 手稿对应位置\n');
fprintf('=============================================================\n');
