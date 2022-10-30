% Copyright (C) 2020-2022 Justin A Fishbone
function plot_save_err_results(test_results, filenameBase, title_txt, xlabel_txt, ylabel_txt)
% plot_save_err_results(test_results, filenameBase, title_txt, xlabel_txt, ylabel_txt)
% Plot and save the simulation results for error vs. k

% Save results
save([filenameBase, '.mat'], 'test_results');

% Plot divergence
genPlot(test_results, filenameBase, title_txt, xlabel_txt, ylabel_txt, test_results.Sq.mean_SP_dkl_fhat_f, test_results.SRocke.mean_SP_dkl_fhat_f, test_results.mmSHR.mean_SP_dkl_fhat_f);
end

function genPlot(test_results, filenameBase, title_txt, xlabel_txt, ylabel_txt, sq_results, sRocke_results, mmSHR_results)

% Constants
FONT_SIZE = 20;
FONT_PARMS = {'FontWeight', 'Bold', 'FontSize', FONT_SIZE};

% Generate plot
figure('Position', [300, 100, 900, 600]);
hold on;
plot(test_results.test_details.k,     sq_results, 'Color', [ 46, 71, 29]/255, 'Marker', '.', 'MarkerSize', 18, 'LineWidth', 2, 'LineStyle', '-');
plot(test_results.test_details.k, sRocke_results, 'Color', [ 47, 85,151]/255, 'Marker', '.', 'MarkerSize', 18, 'LineWidth', 2, 'LineStyle', '-');
plot(test_results.test_details.k,  mmSHR_results, 'Color', [191,144,  0]/255, 'Marker', '.', 'MarkerSize', 18, 'LineWidth', 2, 'LineStyle', '-');
xlabel(xlabel_txt, FONT_PARMS{:});
ylabel(ylabel_txt, FONT_PARMS{:});
grid on;
set(gca, 'Box', 'On', FONT_PARMS{:});
title(title_txt, FONT_PARMS{:});
legend('S-q', 'S-Rocke', 'MM-SHR', 'Location', 'NorthWest');
saveas(gcf, [filenameBase, '.fig']);
saveas(gcf, [filenameBase, '.png']);
end