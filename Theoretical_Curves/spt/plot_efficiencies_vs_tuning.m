% Copyright (C) 2020-2022 Justin A Fishbone
function plot_efficiencies_vs_tuning(qTun, rockeTun, mmTun, yStruct, xLabel, xScale, figName)

figure;
hold on;
plot(linspace(-10,10,length(qTun)), [yStruct.sample], 'LineWidth', 3, 'Color', 'k',             'LineStyle', '--');
plot(qTun, [yStruct.sq],                              'LineWidth', 4, 'Color', [46,71,29]/255);
plot(rockeTun, [yStruct.sRocke],                      'LineWidth', 2, 'Color', [47,85,151]/255);
plot(mmTun, [yStruct.mmShr],                          'LineWidth', 3, 'Color', [191,144,0]/255, 'LineStyle', '-.');
plot(linspace(-10,10,length(qTun)), [yStruct.sBisq],  'LineWidth', 2, 'Color', [197,90,17]/255, 'LineStyle', '--');


xlabel(xLabel);
xlim([-0.4, 1.9]);
set(gca, 'xscale', xScale, 'box', 'on');
ylabel('Asymptotic Relative Efficiency');
ylim([0, 0.85]);
grid on;
legend('SCM', 'S-q', 'S-Rocke', 'MM-SHR', 'S-Bisq');
saveas(gcf, figName);
end