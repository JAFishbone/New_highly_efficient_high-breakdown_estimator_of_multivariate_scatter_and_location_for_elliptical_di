% Copyright (C) 2020-2022 Justin A Fishbone
function plot_save_efficiency_test_results(test_results, filenameBase, title_txt, xlabel_txt, ylabel_txt)
% plot_save_efficiency_test_results(test_results, filenameBase, title_txt, xlabel_txt, ylabel_txt)
% Plot and save the efficiency versus parameter simulation results

% Save results
save([filenameBase, '.mat'], 'test_results');

% Plots
% Gaussian KL divergence efficiency
genPlot(test_results, filenameBase, title_txt, xlabel_txt, ylabel_txt, test_results.Sq.eff_SP_dkl_fhat_f, test_results.SRocke.eff_SP_dkl_fhat_f, test_results.mmSHR.eff_SP_dkl_fhat_f, test_results.Sq2.eff_SP_dkl_fhat_f);
end

function genPlot(test_results, filenameBase, title_txt, xlabel_txt, ylabel_txt, sq_results, sRocke_results, mmSHR_results, sq2_results)

% Constants
FONT_SIZE = 20;
FONT_PARMS = {'FontWeight', 'Bold', 'FontSize', FONT_SIZE};

% Generate plot
figure('Position', [300, 100, 900, 600]);
hold on;
plot(test_results.test_details.distParm,     sq_results, 'Color', [ 46, 71, 29]/255, 'Marker', '.', 'MarkerSize', 18, 'LineWidth', 2, 'LineStyle', '-');
plot(test_results.test_details.distParm, sRocke_results, 'Color', [ 47, 85,151]/255, 'Marker', '.', 'MarkerSize', 18, 'LineWidth', 2, 'LineStyle', '-');
plot(test_results.test_details.distParm,  mmSHR_results, 'Color', [191,144,  0]/255, 'Marker', '.', 'MarkerSize', 18, 'LineWidth', 2, 'LineStyle', '-');
if nargin > 8
    plot(test_results.test_details.distParm,     sq2_results, 'Color', 2.5*[ 46, 71, 29]/255, 'Marker', '.', 'MarkerSize', 18, 'LineWidth', 2, 'LineStyle', '-');
end
ylim([0, 1.05]);
xlabel(xlabel_txt, FONT_PARMS{:});
ylabel(ylabel_txt, FONT_PARMS{:});
grid on;
if min(test_results.test_details.distParm) > 0
    set(gca, 'XScale', 'log');
end
set(gca, 'Box', 'On', FONT_PARMS{:});
title(title_txt, FONT_PARMS{:});
if nargin > 8
    legend('S-q', 'S-Rocke', 'MM-SHR', 'S-q2', 'Location', 'SouthEast');
else
    legend('S-q', 'S-Rocke', 'MM-SHR', 'Location', 'SouthEast');
end
saveas(gcf, [filenameBase, '.fig']);
saveas(gcf, [filenameBase, '.png']);
end