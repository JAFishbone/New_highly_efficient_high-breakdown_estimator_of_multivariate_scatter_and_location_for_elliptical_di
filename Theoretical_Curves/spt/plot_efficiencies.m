% Copyright (C) 2020-2022 Justin A Fishbone
function plot_efficiencies(x, yStruct, xLabel, xScale, figName)

figure;
hold on;
plot(x, [yStruct.sample], 'LineWidth', 3, 'Color', 'k',             'LineStyle', '--');
plot(x, [yStruct.sq],     'LineWidth', 4, 'Color', [46,71,29]/255);
plot(x, [yStruct.sqMax],  'LineWidth', 4, 'Color', [46,71,29]/255, 'LineStyle', ':');
plot(x, [yStruct.sRocke], 'LineWidth', 2, 'Color', [47,85,151]/255);
plot(x, [yStruct.mmShr],  'LineWidth', 3, 'Color', [191,144,0]/255, 'LineStyle', '-.');
plot(x, [yStruct.sBisq],  'LineWidth', 2, 'Color', [197,90,17]/255, 'LineStyle', '--');
plot(x, [yStruct.sqFixed],'LineWidth', 4, 'Color', [25,25,25]/255, 'LineStyle', ':');


xlabel(xLabel);
xlim([x(1),x(end)]);
set(gca, 'xscale', xScale, 'box', 'on');
ylabel('Asymptotic Relative Efficiency');
ylim([0, 1]);
grid on;
legend('SCM', 'S-q', 'S-q Max', 'S-Rocke', 'MM-SHR', 'S-Bisq', 'S-q Fixed');
saveas(gcf, figName);
end