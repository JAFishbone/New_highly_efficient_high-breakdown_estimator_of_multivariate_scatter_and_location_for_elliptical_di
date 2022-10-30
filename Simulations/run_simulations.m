% Copyright (C) 2020-2022 Justin A Fishbone
function run_simulations
% This is a master run control file that calls all of the finite-sample simulations for efficiency and robustness sections

% These all call:
% test_results = estimator_sim(p, n, nsteps, ntrials, distParm, mu, sigma, epsilon, k, rngFun, mleFun, b, sqFun, sqDistParms, q, gam, mmshrCnst, inclBisq, initializeTrue, tol)
% Execute the estimator simulations
%
% INPUTS:
% p              - Dimension
% n              - Number of samples per trial
% nsteps         - Number of steps (length of vector inputs)
% ntrials        - Number of trials per distribution parameter
% distParm       - Distribution parameter value(s) (scalar or vector)
% mu             - px1 true location vector
% sigma          - pxp true scatter matrix
% epsilon        - Contamination proportion
% k              - Contamination value (scalar or vector)
% rngFun         - Function handle to random number generator function, with only input being the distribution parameter (i.e. X = rngFun(distParm(i));) (if there is one)
% mleFun         - Function handle to MLE function, parameterized as: mleFun(X, distParm(i), mu0, sigma0)
% b              - Estimator value b (for maximum breakdown point, set to b=0.5-0.5*p/n)
% sqFun          - S-q estimator tag for sEst() function
% sqDistParms    - Structure(s) for S-q estimator (scalar or vector)
% q              - S-q    tuning parameter  (scalar or vector)
% gam            - Rocke  tuning parameter  (scalar or vector)
% mmshrCnst      - MM-SHR tuning parameter  (scalar or vector)
% inclBisq       - Logical of whether to include S-Bisquare in the run
% initializeTrue - Whether to initialize estimators with true location and scatter
% tol            - Estimator convergence tolerance (Can be empty)
% 


run_fig5b; % Figure 5b -- Max Efficiency Finite-Sample t-Distribution
run_fig8a; % Figure 8a -- Gaussian divergence vs. k, p =  5, n =   5p
run_fig8b; % Figure 8b -- Gaussian divergence vs. k, p = 20, n =   5p
run_fig8c; % Figure 8c -- Gaussian divergence vs. k, p =  5, n = 100p
run_fig8d; % Figure 8d -- Gaussian divergence vs. k, p = 20, n = 100p
end



function run_fig5b % Figure 5b -- Max Efficiency Finite-Sample t-Distribution

% Distribution parameters
NU             = [1, 1.5, 2, 3, 4, 5, 7, 10];
% Simulation Parameters
p              = 20;
nX             = 3;
n              = nX * p;
nsteps         = length(NU);
ntrials        = 6000;
distParm       = NU;
mu             = zeros(p,1);
sigma          = eye(p);
epsilon        = 0;
k              = [];
tol            = 1e-6;
rngFun         = @(nu) PearsonVIIDist.rand(p, n, mu, sigma, nu, (nu+p)/2);
mleFun         = @(X, nu, mu0, sigma0) PearsonVIIDist.mle(X, mu0, sigma0, tol, nu, (nu+p)/2);
b              = 0.5-0.5*p/n; % Maximum breakdown point
sqFun          = 'QT';
sqFixedFun     = 'QCAUCHY';
for i = 1:length(NU)
    sqDistParms(i).nu = NU(i); 
end
q              = 0.999;
gam            = 1;
mmshrCnst      = [1, 1, 1, 1, 1, 1, 1, 1];
inclBisq       = false;
initializeTrue = false;
% Run simulation
results = estimator_sim(p, n, nsteps, ntrials, distParm, mu, sigma, epsilon, k, rngFun, mleFun, b, sqFun, sqDistParms, q, gam, mmshrCnst, inclBisq, initializeTrue, tol, sqFixedFun);

% Plot Results
filenameBase   = 'Max_Eff_vs_Nu-tDist-Finite-p20-n3p';
title_txt      = 'Max Finite-Sample Eff. t-Dist';
xlabel_txt     = 't-Distribution Parameter, \nu';
ylabel_txt     = 'Finite-Sample Relative Efficiency';
plot_save_efficiency_test_results(results, filenameBase, title_txt, xlabel_txt, ylabel_txt);
end



function run_fig8a % Figure 8a -- Gaussian divergence vs. k, p =  5, n =   5p
% Parameters for estimator_sim()
p              = 5;
nX             = 5;
n              = nX * p;
nsteps         = 19;
ntrials        = 250000;
distParm       = NaN; % N/A
mu             = zeros(p,1);
sigma          = eye(p);
epsilon        = 0.1;
k              = [0:0.5:2.5, 2.75:0.25:3.5, 4:0.5:6, 7:10];
tol            = 1e-6;
rngFun         = @(~) GaussDist.rand(p, n, mu, sigma);
mleFun         = @(X, ~, mu0, sigma0) GaussDist.mle(X, mu0, sigma0, tol);
b              = 0.5-0.5*p/n; % Maximum breakdown point
sqFun          = 'QGAUSS';
sqDistParms    = struct;
q              = 0.62;
gam            = 1;
mmshrCnst      = 0.567;
inclBisq       = false;
initializeTrue = false;
% Run simulation
results = estimator_sim(p, n, nsteps, ntrials, distParm, mu, sigma, epsilon, k, rngFun, mleFun, b, sqFun, sqDistParms, q, gam, mmshrCnst, inclBisq, initializeTrue, tol);

% Plot Results
filenameBase   = sprintf('Err_vs_k-Gaussian-p%d-n%d', p, n);
title_txt      = sprintf('Gaussian Error vs. k -- p=%d, n=%d', p, n);
xlabel_txt     = 'k';
ylabel_txt     = 'Small-Sample Divergence';
plot_save_err_results(results, filenameBase, title_txt, xlabel_txt, ylabel_txt);
end


function run_fig8b % Figure 8b -- Gaussian divergence vs. k, p = 20, n =   5p
% Parameters for estimator_sim()
p              = 20;
nX             = 5;
n              = nX * p;
nsteps         = 34;
ntrials        = 20000;
distParm       = NaN; % N/A
mu             = zeros(p,1);
sigma          = eye(p);
epsilon        = 0.1;
k              = [0:0.5:4, 4.2:0.25:8, 8.5:0.5:10, 11:15];
tol            = 1e-6;
rngFun         = @(~) GaussDist.rand(p, n, mu, sigma);
mleFun         = @(X, ~, mu0, sigma0) GaussDist.mle(X, mu0, sigma0, tol);
b              = 0.5-0.5*p/n; % Maximum breakdown point
sqFun          = 'QGAUSS';
sqDistParms    = struct;
q              = 0.919;
gam            = 0.973;
mmshrCnst      = 0.611;
inclBisq       = false;
initializeTrue = false;
% Run simulation
results = estimator_sim(p, n, nsteps, ntrials, distParm, mu, sigma, epsilon, k, rngFun, mleFun, b, sqFun, sqDistParms, q, gam, mmshrCnst, inclBisq, initializeTrue, tol);

% Plot Results
filenameBase   = sprintf('Err_vs_k-Gaussian-p%d-n%d', p, n);
title_txt      = sprintf('Gaussian Error vs. k -- p=%d, n=%d', p, n);
xlabel_txt     = 'k';
ylabel_txt     = 'Small-Sample Divergence';
plot_save_err_results(results, filenameBase, title_txt, xlabel_txt, ylabel_txt);
end


function run_fig8c % Figure 8c -- Gaussian divergence vs. k, p =  5, n = 100p
% Parameters for estimator_sim()
p              = 5;
nX             = 100;
n              = nX * p;
nsteps         = 41;
ntrials        = 4000;
distParm       = NaN; % N/A
mu             = zeros(p,1);
sigma          = eye(p);
epsilon        = 0.1;
k              = 0:0.1:4;
tol            = 1e-6;
rngFun         = @(~) GaussDist.rand(p, n, mu, sigma);
mleFun         = @(X, ~, mu0, sigma0) GaussDist.mle(X, mu0, sigma0, tol);
b              = 0.5-0.5*p/n; % Maximum breakdown point
sqFun          = 'QGAUSS';
sqDistParms    = struct;
q              = 0.537;
gam            = 1;
mmshrCnst      = 0.796;
inclBisq       = false;
initializeTrue = false;
% Run simulation
results = estimator_sim(p, n, nsteps, ntrials, distParm, mu, sigma, epsilon, k, rngFun, mleFun, b, sqFun, sqDistParms, q, gam, mmshrCnst, inclBisq, initializeTrue, tol);

% Plot Results
filenameBase   = sprintf('Err_vs_k-Gaussian-p%d-n%d', p, n);
title_txt      = sprintf('Gaussian Error vs. k -- p=%d, n=%d', p, n);
xlabel_txt     = 'k';
ylabel_txt     = 'Large-Sample Divergence';
plot_save_err_results(results, filenameBase, title_txt, xlabel_txt, ylabel_txt);
end


function run_fig8d % Figure 8d -- Gaussian divergence vs. k, p = 20, n = 100p
% Parameters for estimator_sim()
p              = 20;
nX             = 100;
n              = nX * p;
nsteps         = 71;
ntrials        = 500;
distParm       = NaN; % N/A
mu             = zeros(p,1);
sigma          = eye(p);
epsilon        = 0.1;
k              = 0:0.1:7;
tol            = 1e-6;
rngFun         = @(~) GaussDist.rand(p, n, mu, sigma);
mleFun         = @(X, ~, mu0, sigma0) GaussDist.mle(X, mu0, sigma0, tol);
b              = 0.5-0.5*p/n; % Maximum breakdown point
sqFun          = 'QGAUSS';
sqDistParms    = struct;
q              = 0.909;
gam            = 0.870;
mmshrCnst      = 0.778;
inclBisq       = false;
initializeTrue = false;
% Run simulation
results = estimator_sim(p, n, nsteps, ntrials, distParm, mu, sigma, epsilon, k, rngFun, mleFun, b, sqFun, sqDistParms, q, gam, mmshrCnst, inclBisq, initializeTrue, tol);

% Plot Results
filenameBase   = sprintf('Err_vs_k-Gaussian-p%d-n%d', p, n);
title_txt      = sprintf('Gaussian Error vs. k -- p=%d, n=%d', p, n);
xlabel_txt     = 'k';
ylabel_txt     = 'Large-Sample Divergence';
plot_save_err_results(results, filenameBase, title_txt, xlabel_txt, ylabel_txt);
end

% Ignore uninitialized vector warning:
%#ok<*AGROW>