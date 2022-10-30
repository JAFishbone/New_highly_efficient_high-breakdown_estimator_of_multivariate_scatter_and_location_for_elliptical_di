% Copyright (C) 2020-2022 Justin A Fishbone

% Parameters
P = 5;

% S-q estimator
Q1 = -2;
Q2 = 0.5;
Q3 = 0.9;

% Laplace
LAMBDA = 1;
PSI    = 2;
CHI    = 0;

% Prepare for plotting
close all;
t  = 0:0.001:10;
figure('Position', [300, 100, 900, 600], 'Renderer', 'Painters');
% Colors
C1 = [ 46,  71,  29]/255;
C2 = [ 84, 130,  53]/255;
C3 = [113, 175,  71]/255;

close all;
figure('Position', [0, 100, 1800, 600], 'Renderer', 'Painters');


%% Gaussian Rho
subplot(2,2,1);
hold on;

% Estimator
sq1 = qGauss(Q1);
sq1.pEff = P;
sq2 = qGauss(Q2);
sq2.pEff = P;
sq3 = qGauss(Q3);
sq3.pEff = P;

% Extreme points
[a1g,c1g,rhoA1g,rhoC1g] = sq1.extremesFun();
s1g = rhoC1g - rhoA1g;
[a2g,c2g,rhoA2g,rhoC2g] = sq2.extremesFun();
s2g = rhoC2g - rhoA2g;
[a3g,c3g,rhoA3g,rhoC3g] = sq3.extremesFun();
s3g = rhoC3g - rhoA3g;

% Plot Gaussian rho tilde
plot(t, (sq1.rhoTildFun(t)-rhoA1g)/s1g, 'Color', C1, 'LineWidth', 2.4, 'LineStyle', ':');
plot(t, (sq2.rhoTildFun(t)-rhoA2g)/s2g, 'Color', C2, 'LineWidth', 2.4, 'LineStyle', ':');
plot(t, (sq3.rhoTildFun(t)-rhoA3g)/s3g, 'Color', C3, 'LineWidth', 2.4, 'LineStyle', ':');
% Plot Gaussian rho
plot(t, sq1.rhoFun(t), 'Color', C1, 'LineWidth', 2.5, 'LineStyle', '-');
plot(t, sq2.rhoFun(t), 'Color', C2, 'LineWidth', 2.5, 'LineStyle', '-');
plot(t, sq3.rhoFun(t), 'Color', C3, 'LineWidth', 2.5, 'LineStyle', '-');
% Plot points a and c
scatter([a1g, c1g], [0,1], 600,'Marker', '.', 'MarkerEdgeColor', C1);
scatter([a2g, c2g], [0,1], 600,'Marker', '.', 'MarkerEdgeColor', C2);
scatter([a3g, c3g], [0,1], 600,'Marker', '.', 'MarkerEdgeColor', C3);
% Add text
text(3.65, 0.85, sprintf('$$q=%G$$', Q1), 'Interpreter', 'LaTeX', 'HorizontalAlignment', 'Right', 'FontSize', 22, 'Color', C1);
text(4.10, 0.50, sprintf('$$q=%G$$', Q2), 'Interpreter', 'LaTeX',                                 'FontSize', 22, 'Color', C2);
text(4.35, 0.23, sprintf('$$q=%G$$', Q3), 'Interpreter', 'LaTeX',                                 'FontSize', 22, 'Color', C3);

% Labels
title('Gaussian', 'FontWeight', 'Bold', 'FontSize', 22);
% xlabel('t', 'FontWeight', 'Bold', 'FontSize', 22);
ylabel('Rho', 'FontWeight', 'Bold', 'FontSize', 22);
set(gca, 'FontWeight', 'Bold', 'FontSize', 22, 'Box', 'On');

% Limits
xlim([t(1), t(end)]);
ylim([-0.05, 1.05]);


%% Gaussian W
subplot(2,2,3);
hold on;

% Plot Gaussian S-q weight
plot(t, sq1.wFun(t)/sq1.wFun(P-2), 'Color', C1, 'LineWidth', 2.5, 'LineStyle', '-'); % 0.96*
plot(t, sq2.wFun(t)/sq2.wFun(P-2), 'Color', C2, 'LineWidth', 2.5, 'LineStyle', '-'); % 0.98*
plot(t, sq3.wFun(t)/sq3.wFun(P-2), 'Color', C3, 'LineWidth', 2.5, 'LineStyle', '-'); % 1.02*
% Plot points a and c
scatter([a1g, c1g], [0,0], 600,'Marker', '.', 'MarkerEdgeColor', C1);
scatter([a2g, c2g], [0,0], 600,'Marker', '.', 'MarkerEdgeColor', C2);
scatter([a3g, c3g], [0,0], 600,'Marker', '.', 'MarkerEdgeColor', C3);
% Plot f(d)
plot(t, chi2pdf(t,P)/chi2pdf(P-2,P), 'k', 'LineWidth', 2.5, 'LineStyle', '--');
% Add text
text(4.4, 0.28, sprintf('$$\\leftarrow \\!\\! q=%G$$', Q1),   'Interpreter', 'LaTeX',                                 'FontSize', 22, 'Color', C1);
text(5.8, 0.11, sprintf('$$q=%G \\!\\! \\rightarrow $$', Q2), 'Interpreter', 'LaTeX',                                 'FontSize', 22, 'Color', C2);
text(9.7, 0.87, sprintf('$$q=%G$$', Q3),                      'Interpreter', 'LaTeX', 'HorizontalAlignment', 'Right', 'FontSize', 22, 'Color', C3);
text(9.9, 0.39,         '$$f(t)$$',                           'Interpreter', 'LaTeX', 'HorizontalAlignment', 'Right', 'FontSize', 22);

% Labels
xlabel('\rm\itt', 'FontSize', 22);
ylabel('Weight', 'FontWeight', 'Bold', 'FontSize', 22);
set(gca, 'FontWeight', 'Bold', 'FontSize', 22, 'Box', 'On', 'YTick', []);

% Limits
xlim([t(1), t(end)]);
ylim([-0.05, 1.05]);


%% Laplace Rho
subplot(2,2,2);
hold on;

% Estimator
sq1 = qGenHyper(LAMBDA, PSI, CHI, Q1);
sq1.pEff = P;
sq2 = qGenHyper(LAMBDA, PSI, CHI, Q2);
sq2.pEff = P;
sq3 = qGenHyper(LAMBDA, PSI, CHI, Q3);
sq3.pEff = P;

% Extreme points
[a1l,c1l,rhoA1l,rhoC1l] = sq1.extremesFun();
s1l = rhoC1l - rhoA1l;
[a2l,c2l,rhoA2l,rhoC2l] = sq2.extremesFun();
s2l = rhoC2l - rhoA2l;
[a3l,c3l,rhoA3l,rhoC3l] = sq3.extremesFun();
s3l = rhoC3l - rhoA3l;

% Plot Laplace rho tilde
plot(t, (sq1.rhoTildFun(t)-rhoA1l)/s1l, 'Color', C1, 'LineWidth', 2.4, 'LineStyle', ':');
plot(t, (sq2.rhoTildFun(t)-rhoA2l)/s2l, 'Color', C2, 'LineWidth', 2.4, 'LineStyle', ':');
plot(t, (sq3.rhoTildFun(t)-rhoA3l)/s3l, 'Color', C3, 'LineWidth', 2.4, 'LineStyle', ':');
% Plot Laplace rho
plot(t, sq1.rhoFun(t), 'Color', C1, 'LineWidth', 2.5, 'LineStyle', '-');
plot(t, sq2.rhoFun(t), 'Color', C2, 'LineWidth', 2.5, 'LineStyle', '-');
plot(t, sq3.rhoFun(t), 'Color', C3, 'LineWidth', 2.5, 'LineStyle', '-');
% Plot points a and c
scatter([a1l, c1l], [0,1], 600,'Marker', '.', 'MarkerEdgeColor', C1);
scatter([a2l, c2l], [0,1], 600,'Marker', '.', 'MarkerEdgeColor', C2);
scatter([a3l, c3l], [0,1], 600,'Marker', '.', 'MarkerEdgeColor', C3);
% Add text
text(0.40, 0.90, sprintf('$$q=%G$$', Q1), 'Interpreter', 'LaTeX', 'FontSize', 22, 'Color', C1);
text(1.80, 0.75, sprintf('$$q=%G$$', Q2), 'Interpreter', 'LaTeX', 'FontSize', 22, 'Color', C2);
text(4.00, 0.35, sprintf('$$q=%G$$', Q3), 'Interpreter', 'LaTeX', 'FontSize', 22, 'Color', C3);

% Labels
title('Laplace', 'FontWeight', 'Bold', 'FontSize', 22);
% xlabel('t', 'FontWeight', 'Bold', 'FontSize', 22);
ylabel('Rho', 'FontWeight', 'Bold', 'FontSize', 22);
set(gca, 'FontWeight', 'Bold', 'FontSize', 22, 'Box', 'On');

% Limits
xlim([t(1), t(end)]);
ylim([-0.05, 1.05]);


%% Laplace Weights
subplot(2,2,4);
hold on;


% Plot Laplace S-q weight
plot(t, sq1.wFun(t)/sq1.wFun(t(2)), 'Color', C1, 'LineWidth', 2.5, 'LineStyle', '-');
plot(t, sq2.wFun(t)/sq2.wFun(t(2)), 'Color', C2, 'LineWidth', 2.5, 'LineStyle', '-');
plot(t, sq3.wFun(t)/sq3.wFun(t(2)), 'Color', C3, 'LineWidth', 2.5, 'LineStyle', '-');
% Plot points a and c
scatter([a1l, c1l], [0,0], 600,'Marker', '.', 'MarkerEdgeColor', C1);
scatter([a2l, c2l], [0,0], 600,'Marker', '.', 'MarkerEdgeColor', C2);
scatter([a3l, c3l], [0,0], 600,'Marker', '.', 'MarkerEdgeColor', C3);
% Plot Density
plot(t, t.^(P/2-1).*sqrt(2*t).^(1-P/2).*besselk(1-P/2,sqrt(2*t)) / (t(2)^(P/2-1).*sqrt(2*t(2)).^(1-P/2).*besselk(1-P/2,sqrt(2*t(2)))) , 'k', 'LineWidth', 2.5, 'LineStyle', '--');
% Add text
% text(0.25, 0.85, sprintf('$$q=%G$$', Q1),                    'Interpreter', 'LaTeX', 'FontSize', 22, 'Color', C1);
% text(1.80, 0.75, sprintf('$$q=%G$$', Q2),                    'Interpreter', 'LaTeX', 'FontSize', 22, 'Color', C2);
text(0.35, 0.50, sprintf('$$\\leftarrow \\!\\! q=%G$$', Q3), 'Interpreter', 'LaTeX', 'FontSize', 22, 'Color', C3);
text(0.40, 0.90,         '$$f(t)$$',                         'Interpreter', 'LaTeX', 'FontSize', 22);

% Labels
xlabel('\rm\itt', 'FontSize', 22);
ylabel('Weight', 'FontWeight', 'Bold', 'FontSize', 22);
set(gca, 'FontWeight', 'Bold', 'FontSize', 22, 'Box', 'On', 'YTick', []);

% Limits
xlim([t(1), t(end)]);
ylim([-0.05, 1.05]);

%%
saveas(gcf, 'rho_and_w_functions.png');
print('-painters', 'rho_and_w_functions.eps', '-depsc');

