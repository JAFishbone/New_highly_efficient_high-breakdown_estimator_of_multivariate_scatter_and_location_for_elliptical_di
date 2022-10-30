% Copyright (C) 2020-2022 Justin A Fishbone

% Parameters
P        = 10;
% Efficiency tuning
Q        = 0.87;
SIG_Q    = 0.842;
ROCKE_G  = 1;
% Consistency parameter, sigma
SIG_MM   = 2.55; % c*S
SIG_ROCK = 9.75;
SIG_BISQ = 46;

% Plotting Parameters
d  = 0:0.01:25;
FONT_SIZE = 20;
FONT_PARMS = {'FontWeight', 'Bold', 'FontSize', FONT_SIZE};

close all;
figure('Position', [300, 50, 900, 600], 'Renderer', 'Painters');
hold on;

% Generate the curves to plot
% PDF
pdf_d   = chi2pdf(d, P)/chi2pdf(P-2, P);
% MLE
w_mle   = ones(size(d));
% Bisquare
w_bisq  = (1-(d/SIG_BISQ)).^2 .* (d/SIG_BISQ <= 1);
% Rocke
w_rocke = (1 - ((d/SIG_ROCK-1)/ROCKE_G).^2) .* ( 1-ROCKE_G <= d/SIG_ROCK & d/SIG_ROCK <= 1+ROCKE_G);
% MM-SHR
y  = d/SIG_MM;
yy = -1.944+1.728*y-0.312*y.^2+0.016*y.^3;
w_mmshr        = yy;
w_mmshr(y>9)  = 0;
w_mmshr(y<=4) = 1;
% S-q
sq = qGauss(Q);
sq.pEff = P;
w_q = sq.wFun(d/SIG_Q)/max(sq.wFun(d/SIG_Q));

% Colors
cSq = [ 46,  71,  29]/255;
cSr = [ 47,  85, 151]/255;
cMM = [191, 144,   0]/255;
cSb = [197,  90,  17]/255; 

% Plot curves
plot(d, w_mle,   'Color', 'k', 'LineWidth', 3, 'LineStyle', '-');
plot(d, pdf_d,   'Color', 'k', 'LineWidth', 4, 'LineStyle', ':');
plot(d, w_bisq,  'Color', cSb, 'LineWidth', 2, 'LineStyle', '--');
plot(d, w_mmshr, 'Color', cMM, 'LineWidth', 5, 'LineStyle', '-.');
plot(d, w_rocke, 'Color', cSr, 'LineWidth', 2, 'LineStyle', '-');
plot(d, w_q,     'Color', cSq, 'LineWidth', 4, 'LineStyle', '-');

% Add text
text(23.00, 0.97, '$$w_\mathrm{mle}$$',                     'Interpreter', 'LaTeX',                                 'FontSize', FONT_SIZE);
text(03.90, 0.40, '$$\leftarrow \!\! f(d)$$',               'Interpreter', 'LaTeX',                                 'FontSize', FONT_SIZE);
text(22.00, 0.29, '$$\leftarrow \! w_\mathrm{bisq}$$',      'Interpreter', 'LaTeX',                                 'FontSize', FONT_SIZE, 'Color', cSb);
text(13.35, 0.82, '$$\rightarrow$$', 'Interpreter', 'LaTeX', 'HorizontalAlignment', 'Right', 'FontSize', FONT_SIZE, 'Color', cMM);
text(13.40, 0.79, '$$w_\mathrm{\scriptscriptstyle{SHR}}$$', 'Interpreter', 'LaTeX', 'HorizontalAlignment', 'Right', 'FontSize', FONT_SIZE, 'Color', cMM);
text(14.91, 0.74, '$$\leftarrow \! w_{\gamma}$$',           'Interpreter', 'LaTeX', 'HorizontalAlignment', 'Left',  'FontSize', FONT_SIZE, 'Color', cSr);
text(21.50, 0.15, '$$\leftarrow \! w_q$$',                  'Interpreter', 'LaTeX',                                 'FontSize', FONT_SIZE, 'Color', cSq);

% Limits
xlim([d(1), d(end)]);
ylim([-0.05, 1.05]);

% Labels
xlabel('Squared Distance, \rm\itd', FONT_PARMS{:});
ylabel('Weight', FONT_PARMS{:});
set(gca, 'Box', 'On', FONT_PARMS{:});

%%
saveas(gcf, 'w_functions.png');
print('-painters', 'w_functions.eps', '-depsc');