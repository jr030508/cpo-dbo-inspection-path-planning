%% ============================================================
%  图6 重绘：30节点基础场景 四算法收敛曲线 (a)DBO (b)IPSO (c)CPO-DBO (d)CPO
%  数据来源：从论文嵌入原图逐像素提取，平台值锚定 51.4888 a.u.
%  输出：output/fig6_convergence_cost.png (300dpi)
%% ============================================================
scriptDir = fileparts(mfilename('fullpath'));
cd(scriptDir);
outDir = fullfile(scriptDir, 'output');
if ~exist(outDir, 'dir'); mkdir(outDir); end
FN = 'Times New Roman';

S = load(fullfile(scriptDir, 'data', 'fig6_curves.mat'));
cA = S.conv_a;  cB = S.conv_b;  cC = S.conv_c;  cD = S.conv_d;
x  = (1:5000)';
dtX = S.dtX;   % datatip X 位置 [444;939;159;312]

titles = {'(a) DBO', '(b) IPSO', '(c) CPO-DBO', '(d) CPO'};
curves = {cA, cB, cC, cD};
ylims  = {[50 110], [50 120], [50 110], [50 110]};
dtCol  = {444, 939, 159, 312};

f = figure('Color','w','Position',[60,60,1080,560],'Visible','off');
for p = 1:4
    subplot(1,4,p); hold on; box on;
    v = curves{p};
    plot(x, v, '-', 'LineWidth', 1.4, 'Color', [0 0.4470 0.7410]);
    xlim([0 5000]); ylim(ylims{p});
    set(gca, 'XTick', 0:500:5000, 'YTick', 50:10:ylims{p}(2), ...
        'FontName', FN, 'FontSize', 9, 'LineWidth', 0.9);
    xlabel('Iterations', 'FontName', FN, 'FontSize', 11);
    ylabel('Optimal exposure cost (a.u.)', 'FontName', FN, 'FontSize', 11);
    % 图例（右下小字）
    legend('Optimal exposure cost', 'FontName', FN, 'FontSize', 8, ...
        'Location', 'northeast');
    % datatip 复刻：蓝色文字框 + 标记点
    xv = dtCol{p};
    yv = v(max(1,min(5000,round(xv))));
    plot(xv, yv, '.', 'MarkerSize', 9, 'Color', [0 0.3 0.6]);
    annotation('textbox', [0.125+ (p-1)*0.245 + 0.10, 0.30, 0.10, 0.10], ...
        'String', sprintf('X: %d\nY: %.4f', xv, yv), ...
        'FitBoxToText','on', 'BackgroundColor','w', 'EdgeColor', [0.6 0.6 0.6], ...
        'Color', [0 0.3 0.6], 'FontName', FN, 'FontSize', 7.5);
    % 面板字母
    text('Units','normalized', 'Position',[-0.28 1.03 0], 'String', titles{p}, ...
        'FontName', FN, 'FontSize', 13, 'FontWeight', 'bold');
end

print(f, '-dpng', '-r300', fullfile(outDir, 'fig6_convergence_cost.png'));
close(f);
fprintf('[图6] 完成 -> output/fig6_convergence_cost.png\n');
