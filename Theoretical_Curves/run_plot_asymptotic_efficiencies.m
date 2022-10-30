% Copyright (C) 2020-2022 Justin A Fishbone

% Run parameters
q = 0.998;
gam = 1;
mmPrm = [];
p = 20;
warning('off');
% Numerically differentiate for simplicity
DEL = 1e-7;
dFun = @(x, fun) (fun(x+DEL) - fun(x-DEL)) ./ (2*DEL);
tic



%% Run t-distribution plot
NU = 0.9:0.1:10;
maxInt = inf;
clear efficienciesT;
for i = 1:length(NU)
    nu = NU(i);
    beta = gamma((nu+p)/2)/gamma((nu+p)/2-p/2)/nu^(p/2)/gamma(p/2);
    fd = @(d) beta * d.^(p/2-1) .* (1+d/nu).^(-(nu+p)/2);
    phi = @(d) (1+d/nu).^(-(nu+p)/2);
    dphi = @(d) -(nu+p)/(2*nu) * (1+d/nu).^(-(nu+p)/2-1);
 
    % Estimators
    sqFun = qPearsonVII(nu,(nu+p)/2);
    sqFunFixed = qPearsonVII(1,(1+p)/2);
    wMLE = @(d) (nu+p)./(nu+d);
    
    efficienciesT(i) = calc_efficiencies(p, beta, fd, dphi, sqFun, sqFunFixed, wMLE, maxInt, q, gam, mmPrm); %#ok<SAGROW>
end
plot_efficiencies(NU, efficienciesT, 't-Distribution Parameter, \nu', 'log', 'MaxEff_vs_tDist.fig');



%% Run Kotz type distribution plot
N = 0; r = 0.5; 
S = 2.^(-1.75:0.05:1.60);
maxInt = inf;
clear efficienciesKotz;
for i = 1:length(S)
    s = S(i);

    beta = s/gamma((2*N+p)/2/s)*r^((2*N+p)/2/s);
    fd = @(d) beta * d.^(p/2-1) .* d.^N .* exp(-r*d.^s);
    phi = @(d) d.^N .* exp(-r*d.^s);
    dphi = @(d) N*d.^(N-1) .* exp(-r*d.^(s)) - r*s*d.^(N+s-1) .* exp(-r*d.^s);

    % Estimators
    sqFun = qKotz(s,N,r);
    sqFunFixed = qKotz(1,N,r);
    wMLE = @(d) -2 * ( N*d.^(-1) - r*s*d.^(s-1) );
    
    efficienciesKotz(i) = calc_efficiencies(p, beta, fd, dphi, sqFun, sqFunFixed, wMLE, maxInt, q, gam, mmPrm); %#ok<SAGROW>
end
plot_efficiencies(S, efficienciesKotz, 'Kotz Type Parameter, s', 'log', 'MaxEff_vs_KotzType.fig');



%% Run Variance Gamma Plot
chi = 0; psi = 2;
LAMBDA = 2.^(-1.4:0.1:4.4);
maxInt = 1e5;
clear efficienciesVG;
for i = 1:length(LAMBDA)
    lambda = LAMBDA(i);
    if chi==0
        beta = psi^(p/4+lambda/2)/2^(p/2+lambda-1)/gamma(lambda)/gamma(p/2);
    else
        beta = psi^(p/4)/2^(p/2)/chi^(lambda/2)/besselk(lambda,sqrt(chi*psi))/gamma(p/2);
    end
    fd = @(d) beta * d.^(p/2-1) .* besselk(lambda-p/2, sqrt(psi*(chi+d)))./(chi+d).^(p/4-lambda/2);
    phi = @(d) besselk(lambda-p/2, sqrt(psi*(chi+d)))./(chi+d).^(p/4-lambda/2);
    dphi = @(d) dFun(d, phi);
 
    % Estimators
    sqFun = qGenHyper(lambda,psi,chi);
    sqFunFixed = qGenHyper(1,psi,chi);
    wMLE = @(d) sqrt(psi./(chi+d)).*besselk(lambda-p/2+1, sqrt(psi*(chi+d)))./besselk(lambda-p/2, sqrt(psi*(chi+d))) - 2*(lambda-p/2)./(chi+d); %dphi(d) ./ phi(d);
    
    efficienciesVG(i) = calc_efficiencies(p, beta, fd, dphi, sqFun, sqFunFixed, wMLE, maxInt, q, gam, mmPrm); %#ok<SAGROW>
end
plot_efficiencies(LAMBDA, efficienciesVG, 'Variance Gamma Parameter, \lambda', 'log', 'MaxEff_vs_VarianceGamma.fig');



%% Run Generalized Hyperbolic Plot
chi = 1; psi = 2;
LAMBDA = -15:0.1:15;
maxInt = 1e5;
clear efficienciesGH;
for i = 1:length(LAMBDA)
    lambda = LAMBDA(i);
    if chi==0
        beta = psi^(p/4+lambda/2)/2^(p/2+lambda-1)/gamma(lambda)/gamma(p/2);
    else
        beta = psi^(p/4)/2^(p/2)/chi^(lambda/2)/besselk(lambda,sqrt(chi*psi))/gamma(p/2);
    end
    fd = @(d) beta * d.^(p/2-1) .* besselk(lambda-p/2, sqrt(psi*(chi+d)))./(chi+d).^(p/4-lambda/2);
    phi = @(d) besselk(lambda-p/2, sqrt(psi*(chi+d)))./(chi+d).^(p/4-lambda/2);
    dphi = @(d) dFun(d, phi);
 
    % Estimators
    sqFun = qGenHyper(lambda,psi,chi);
    sqFunFixed = qGenHyper(-0.5,psi,chi);
    wMLE = @(d) sqrt(psi./(chi+d)).*besselk(lambda-p/2+1, sqrt(psi*(chi+d)))./besselk(lambda-p/2, sqrt(psi*(chi+d))) - 2*(lambda-p/2)./(chi+d); %dphi(d) ./ phi(d);
    
    efficienciesGH(i) = calc_efficiencies(p, beta, fd, dphi, sqFun, sqFunFixed, wMLE, maxInt, q, gam, mmPrm); %#ok<SAGROW>
end
plot_efficiencies(LAMBDA, efficienciesGH, 'Generalized Hyperbolic Parameter, \lambda', 'linear', 'MaxEff_vs_GeneralizedHyperbolic.fig');




%% Run t-distribution asymptotic efficiency vs. tuning
NSTEPS = 100;
qTun     = linspace(-0.5, 0.999, NSTEPS);
rockeTun = linspace( 0.01, 1,     NSTEPS);
mmTun    = linspace( 0.01, 2,     NSTEPS);
bisqTun  = linspace(   0, 0,     NSTEPS);
maxInt = inf;
clear efficienciesT_vs_tun;
for i = 1:NSTEPS
    % Distribution
    nu = 1;
    beta = gamma((nu+p)/2)/gamma((nu+p)/2-p/2)/nu^(p/2)/gamma(p/2);
    fd = @(d) beta * d.^(p/2-1) .* (1+d/nu).^(-(nu+p)/2);
    phi = @(d) (1+d/nu).^(-(nu+p)/2);
    dphi = @(d) -(nu+p)/(2*nu) * (1+d/nu).^(-(nu+p)/2-1);
 
    % Tuning parameters
    q = qTun(i);
    gam = rockeTun(i);
    mmPrm = mmTun(i);
        
    % Estimators
    sqFun = qPearsonVII(nu,(nu+p)/2);
    sqFunFixed = qPearsonVII(1,(1+p)/2);
    wMLE = @(d) (nu+p)./(nu+d);
    
    efficienciesT_vs_tun(i) = calc_efficiencies(p, beta, fd, dphi, sqFun, sqFunFixed, wMLE, maxInt, q, gam, mmPrm); %#ok<SAGROW>
end
plot_efficiencies_vs_tuning(qTun, rockeTun, mmTun, efficienciesT_vs_tun, 'Tuning Parameter', 'linear', 'tDist_efficiency_vs_tun.fig');


%%
toc
warning('on');