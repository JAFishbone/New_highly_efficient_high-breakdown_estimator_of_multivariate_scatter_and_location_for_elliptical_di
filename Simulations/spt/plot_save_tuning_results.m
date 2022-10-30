% Copyright (C) 2020-2022 Justin A Fishbone
function plot_save_tuning_results(test_results, filenameBase, title_txt, xlabel_txt, ylabel_txt)
% plot_save_tuning_results(test_results, filenameBase, title_txt, xlabel_txt, ylabel_txt)
% Plot and save the simulation results for efficiency vs. tuning

% Save results
save([filenameBase, '.mat'], 'test_results');

% Plots
% Gaussian squared Hellinger distance efficiency
genPlot(test_results, [filenameBase, '-H2'],  [title_txt, ' - H2 Efficiency'],  xlabel_txt, ylabel_txt, test_results.Sq.eff_SP_h2        , test_results.SRocke.eff_SP_h2        , test_results.mmSHR.eff_SP_h2        , test_results.SBisq.eff_SP_h2);
% Gaussian KL divergence efficiency
genPlot(test_results, [filenameBase, '-DKL'], [title_txt, ' - DKL Efficiency'], xlabel_txt, ylabel_txt, test_results.Sq.eff_SP_dkl_fhat_f, test_results.SRocke.eff_SP_dkl_fhat_f, test_results.mmSHR.eff_SP_dkl_fhat_f, test_results.SBisq.eff_SP_dkl_fhat_f);
% xi efficiency
genPlot(test_results, [filenameBase, '-xi'],  [title_txt, ' - xi Efficiency'],  xlabel_txt, ylabel_txt, test_results.Sq.eff_SP_xi        , test_results.SRocke.eff_SP_xi        , test_results.mmSHR.eff_SP_xi        , test_results.SBisq.eff_SP_xi);
end

function genPlot(test_results, filenameBase, title_txt, xlabel_txt, ylabel_txt, sq_results, sRocke_results, mmSHR_results, sBisq_results)

% Constants
FONT_SIZE = 20;
FONT_PARMS = {'FontWeight', 'Bold', 'FontSize', FONT_SIZE};

% Generate plot
figure('Position', [300, 100, 900, 600]);
hold on;
plot(test_results.test_details.q,                sq_results, 'Color', [ 46, 71, 29]/255, 'Marker', '.', 'MarkerSize', 18, 'LineWidth', 2, 'LineStyle', '-');
plot(test_results.test_details.gam,          sRocke_results, 'Color', [ 47, 85,151]/255, 'Marker', '.', 'MarkerSize', 18, 'LineWidth', 2, 'LineStyle', '-');
plot(test_results.test_details.mmshrCnst,     mmSHR_results, 'Color', [191,144,  0]/255, 'Marker', '.', 'MarkerSize', 18, 'LineWidth', 2, 'LineStyle', '-');
plot([-100,100], [mean(sBisq_results), mean(sBisq_results)], 'Color', [197, 90, 17]/255,                                  'LineWidth', 2, 'LineStyle', '--');
ylim([0, 0.85]);
xlim([min(test_results.test_details.q)-0.1, max(test_results.test_details.mmshrCnst)+0.1]);
xlabel(xlabel_txt, FONT_PARMS{:});
ylabel(ylabel_txt, FONT_PARMS{:});
grid on;
set(gca, 'Box', 'On', FONT_PARMS{:});
title(title_txt, FONT_PARMS{:});
legend('S-q', 'S-Rocke', 'MM-SHR', 'S-Bisq', 'Location', 'SouthEast');
saveas(gcf, [filenameBase, '.fig']);
saveas(gcf, [filenameBase, '.png']);
end