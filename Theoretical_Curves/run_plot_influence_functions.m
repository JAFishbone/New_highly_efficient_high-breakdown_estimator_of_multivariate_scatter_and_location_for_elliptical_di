% Copyright (C) 2020-2022 Justin A Fishbone

%%% p=10 Gaussian with MM & S-q at 80%, S-Rocke at 77% efficiency
p = 10;
phi = @(d) exp(-d/2);
dphi = @(d) -0.5*exp(-d/2);
beta = 1/gamma(p/2)*0.5^(p/2);
fd = @(d) beta * d.^(p/2-1) .* exp(-0.5*d);


dz = 0:0.01:25;

% S-q
q = 0.87;
sFun = qGauss;
sFun.pEff = p;  sFun.tun = q;
%
rhoFun = @(t) sFun.rhoFun(t);
wFun = @(t) sFun.wFun(t);
sigma = fzero( @(sigma) integral( @(d) fd(d).*rhoFun(d/sigma), 0, inf ) - 0.5, [1e-6, 1e3]);
lambda1 = -2*beta* integral( @(d) 1/sigma * d.^(p/2+1) .* wFun(d/sigma) .* dphi(d), 0, inf );
%
alphaSQ = p*(p+2)*(dz/sigma).*sFun.wFun(dz/sigma)/lambda1;

% S-Rocke
gam = 1;
sFun = sRocke;
sFun.pEff = p;  sFun.tun = gam;
%
rhoFun = @(t) sFun.rhoFun(t);
wFun = @(t) sFun.wFun(t);
sigma = fzero( @(sigma) integral( @(d) fd(d).*rhoFun(d/sigma), 0, inf ) - 0.5, [1e-6, 1e3]);
lambda1 = -2*beta* integral( @(d) 1/sigma * d.^(p/2+1) .* wFun(d/sigma) .* dphi(d), 0, inf );
%
alphaRCK = p*(p+2)*(dz/sigma).*sFun.wFun(dz/sigma)/lambda1;


% S-Bisq
sFun = sBisq;
sFun.pEff = p;
%
rhoFun = @(t) sFun.rhoFun(t);
wFun = @(t) sFun.wFun(t);
sigma = fzero( @(sigma) integral( @(d) fd(d).*rhoFun(d/sigma), 0, inf ) - 0.5, [1e-6, 1e3]);
lambda1 = -2*beta* integral( @(d) 1/sigma * d.^(p/2+1) .* wFun(d/sigma) .* dphi(d), 0, inf );
%
alphaBisq = p*(p+2)*(dz/sigma).*sFun.wFun(dz/sigma)/lambda1;


% MM-SHR
sigmaMM = 2.55;
%
wMM = @(t) ( (-1.944+1.728*t-0.312*t.^2+0.016*t.^3) .* (t<=9 & t>4) + (t<=4) );
dwMM = @(t) ( (1.728-0.624*t+0.048*t.^2) .* (t<=9 & t>4) );
lambda1 = -2*beta* integral( @(d) 1/sigmaMM * d.^(p/2+1) .* wMM(d/sigmaMM) .* dphi(d), 0, inf );
alphaMM = p*(p+2)*(dz/sigmaMM).*wMM(dz/sigmaMM)/lambda1;

  

% SCM
lambda1  = integral( @(d) (p+2)*fd(d).*d, 0, inf);
alphaSCM = p*(p+2)*dz/lambda1;



%% Plot results
figure;
hold on;


plot(dz, alphaSCM, 'LineWidth', 3, 'Color', 'k',             'LineStyle', '--');
plot(dz, alphaSQ,  'LineWidth', 4, 'Color', [46,71,29]/255);
plot(dz, alphaRCK, 'LineWidth', 2, 'Color', [47,85,151]/255);
plot(dz, alphaMM,  'LineWidth', 3, 'Color', [191,144,0]/255, 'LineStyle', '-.');
plot(dz, alphaBisq,'LineWidth', 2, 'Color', [197,90,17]/255, 'LineStyle', '--');


xlabel('d_z');
xlim([dz(1), dz(end)]);
set(gca, 'box', 'on');
ylabel('\alpha_{\bf\Sigma}(d_z) from Influence Function');
ylim([0, 35/2]);
grid on;
legend('SCM', 'S-q', 'S-Rocke', 'MM-SHR', 'S-Bisq');
saveas(gcf, 'influence_functions.fig');
