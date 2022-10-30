% Copyright (C) 2020-2022 Justin A Fishbone
function run_computational_sims
% Run the stability and computational efficiency simulations from Section 5

% These call:
% test_results = stability_sim(p, n, ntrials, epsilon, k, rngFun, mleFun, b, sqFun, sqDistParms, q, gam, mmshrCnst, tol, initFrac)
% test_results = computational_efficiency_sim(p, n, ntrials, epsilon, k, rngFun, b, sqFun, sqDistParms, q, gam, mmshrCnst, tol)
%

%
% test_results = stability_sim(p, n, ntrials, epsilon, k, rngFun, mleFun, b, sqFun, sqDistParms, q, gam, mmshrCnst, tol, initFrac)
% Execute the estimator stability simulations
%
% INPUTS:
% p              - Dimension
% n              - Number of samples per trial
% ntrials        - Number of trials per distribution parameter
% epsilon        - Contamination proportion
% k              - Contamination value (vector length 3 corresponding to each estimator: [S-q, S-Rocke, MM-SHR])
% rngFun         - Function handle to random number generator function
% mleFun         - Function handle to MLE function, parameterized as: mleFun(X)
% b              - Estimator value b (for maximum breakdown point, set to b=0.5-0.5*p/n)
% sqFun          - S-q estimator tag for sEst() function
% sqDistParms    - Structure for S-q estimator
% q              - S-q    tuning parameter
% gam            - Rocke  tuning parameter
% mmshrCnst      - MM-SHR tuning parameter
% tol            - Estimator convergence tolerance (Can be empty)
% initFrac       - Proportion of n samples that are used for MLE initialization of estimators
% 

% test_results = computational_efficiency_sim(p, n, ntrials, epsilon, k, rngFun, b, sqFun, sqDistParms, q, gam, mmshrCnst, tol)
% Execute the estimator stability simulations
%
% INPUTS:
% p              - Dimension
% n              - Number of samples per trial
% ntrials        - Number of trials per distribution parameter
% epsilon        - Contamination proportion
% k              - Contamination value (vector length 3 corresponding to each estimator: [S-q, S-Rocke, MM-SHR])
% rngFun         - Function handle to random number generator function
% b              - Estimator value b (for maximum breakdown point, set to b=0.5-0.5*p/n)
% sqFun          - S-q estimator tag for sEst() function
% sqDistParms    - Structure for S-q estimator
% q              - S-q    tuning parameter
% gam            - Rocke  tuning parameter
% mmshrCnst      - MM-SHR tuning parameter
% tol            - Estimator convergence tolerance (Can be empty)
% 

% Constants across all runs
result_file = 'computational_simulation_results.txt';
mleFun      = @(X) GaussDist.mle(X);
sqFun       = 'QGAUSS';
sqDistParms = struct();
tol         = 1e-10;
initFrac    = 0.25;
muFun       = @(p) zeros(p,1);
sigmaFun    = @(p) eye(p);

% Prepare for the runs
stability_results = {};
comp_eff_results = {};
t0 = now;
% For reproducibility, seed the random number generator
rng('default');


%% P=5, n=5p, EPSILON=0%
p         = 5;
n         = 5*p;
ntrials   = 1000;
epsilon   = 0;
k         = [3.0, 3.25, 2.75]; % [MLE, S-q, S-Rocke, MM-SHR] - S-q, S-Rocke, MM-SHR are set to max from error sim
rngFun    = @(~) GaussDist.rand(p, n, muFun(p), sigmaFun(p));
b         = 0.5-0.5*p/n; % Maximum breakdown point
q         = 0.62;
gam       = 1;
mmshrCnst = 0.567;
stability_results{end+1} = stability_sim(               p, n, ntrials, epsilon, k, rngFun, mleFun, b, sqFun, sqDistParms, q, gam, mmshrCnst, tol, initFrac);
comp_eff_results{end+1}  = computational_efficiency_sim(p, n, ntrials, epsilon, k, rngFun,         b, sqFun, sqDistParms, q, gam, mmshrCnst, tol);

%% P=5, n=100p, EPSILON=0%
p         = 5;
n         = 100*p;
ntrials   = 1000;
epsilon   = 0;
k         = [2.7, 2.8, 2.75]; % [S-q, S-Rocke, MM-SHR] - S-q, S-Rocke, MM-SHR are set to max from error sim
rngFun    = @(~) GaussDist.rand(p, n, muFun(p), sigmaFun(p));
b         = 0.5-0.5*p/n; % Maximum breakdown point
q         = 0.537;
gam       = 1;
mmshrCnst = 0.796;
stability_results{end+1} = stability_sim(               p, n, ntrials, epsilon, k, rngFun, mleFun, b, sqFun, sqDistParms, q, gam, mmshrCnst, tol, initFrac);
comp_eff_results{end+1}  = computational_efficiency_sim(p, n, ntrials, epsilon, k, rngFun,         b, sqFun, sqDistParms, q, gam, mmshrCnst, tol);


%% P=20, n=5p, EPSILON=0%
p         = 20;
n         = 5*p;
ntrials   = 1000;
epsilon   = 0;
k         = [6.1, 6.0, 6.0]; % [S-q, S-Rocke, MM-SHR] - S-q, S-Rocke, MM-SHR are set to max from error sim
rngFun    = @(~) GaussDist.rand(p, n, muFun(p), sigmaFun(p));
b         = 0.5-0.5*p/n; % Maximum breakdown point
q         = 0.919;
gam       = 0.973;
mmshrCnst = 0.611;
stability_results{end+1} = stability_sim(               p, n, ntrials, epsilon, k, rngFun, mleFun, b, sqFun, sqDistParms, q, gam, mmshrCnst, tol, initFrac);
comp_eff_results{end+1}  = computational_efficiency_sim(p, n, ntrials, epsilon, k, rngFun,         b, sqFun, sqDistParms, q, gam, mmshrCnst, tol);


%% P=20, n=100p, EPSILON=0%
p         = 20;
n         = 100*p;
ntrials   = 1000;
epsilon   = 0;
k         = [5.8, 5.3, 5.7]; % [S-q, S-Rocke, MM-SHR] - S-q, S-Rocke, MM-SHR are set to max from error sim
rngFun    = @(~) GaussDist.rand(p, n, muFun(p), sigmaFun(p));
b         = 0.5-0.5*p/n; % Maximum breakdown point
q         = 0.909;
gam       = 0.870;
mmshrCnst = 0.778;
stability_results{end+1} = stability_sim(               p, n, ntrials, epsilon, k, rngFun, mleFun, b, sqFun, sqDistParms, q, gam, mmshrCnst, tol, initFrac);
comp_eff_results{end+1}  = computational_efficiency_sim(p, n, ntrials, epsilon, k, rngFun,         b, sqFun, sqDistParms, q, gam, mmshrCnst, tol);


%% P=5, n=5p, EPSILON=10%
p         = 5;
n         = 5*p;
ntrials   = 1000;
epsilon   = 0.1;
k         = [3.0, 3.25, 2.75]; % [S-q, S-Rocke, MM-SHR] - S-q, S-Rocke, MM-SHR are set to max from error sim
rngFun    = @(~) GaussDist.rand(p, n, muFun(p), sigmaFun(p));
b         = 0.5-0.5*p/n; % Maximum breakdown point
q         = 0.620;
gam       = 1;
mmshrCnst = 0.567;
stability_results{end+1} = stability_sim(               p, n, ntrials, epsilon, k, rngFun, mleFun, b, sqFun, sqDistParms, q, gam, mmshrCnst, tol, initFrac);
comp_eff_results{end+1}  = computational_efficiency_sim(p, n, ntrials, epsilon, k, rngFun,         b, sqFun, sqDistParms, q, gam, mmshrCnst, tol);


%% P=5, n=100p, EPSILON=10%
p         = 5;
n         = 100*p;
ntrials   = 1000;
epsilon   = 0.1;
k         = [2.7, 2.8, 2.75]; % [S-q, S-Rocke, MM-SHR] - S-q, S-Rocke, MM-SHR are set to max from error sim
rngFun    = @(~) GaussDist.rand(p, n, muFun(p), sigmaFun(p));
b         = 0.5-0.5*p/n; % Maximum breakdown point
q         = 0.537;
gam       = 1;
mmshrCnst = 0.796;
stability_results{end+1} = stability_sim(               p, n, ntrials, epsilon, k, rngFun, mleFun, b, sqFun, sqDistParms, q, gam, mmshrCnst, tol, initFrac);
comp_eff_results{end+1}  = computational_efficiency_sim(p, n, ntrials, epsilon, k, rngFun,         b, sqFun, sqDistParms, q, gam, mmshrCnst, tol);


%% P=20, n=5p, EPSILON=10%
p         = 20;
n         = 5*p;
ntrials   = 1000;
epsilon   = 0.1;
k         = [6.1, 6.0, 6.0]; % [S-q, S-Rocke, MM-SHR] - S-q, S-Rocke, MM-SHR are set to max from error sim
rngFun    = @(~) GaussDist.rand(p, n, muFun(p), sigmaFun(p));
b         = 0.5-0.5*p/n; % Maximum breakdown point
q         = 0.919;
gam       = 0.973;
mmshrCnst = 0.611;
stability_results{end+1} = stability_sim(               p, n, ntrials, epsilon, k, rngFun, mleFun, b, sqFun, sqDistParms, q, gam, mmshrCnst, tol, initFrac);
comp_eff_results{end+1}  = computational_efficiency_sim(p, n, ntrials, epsilon, k, rngFun,         b, sqFun, sqDistParms, q, gam, mmshrCnst, tol);


%% P=20, n=100p, EPSILON=10%
p         = 20;
n         = 100*p;
ntrials   = 1000;
epsilon   = 0.1;
k         = [5.8, 5.3, 5.7]; % [S-q, S-Rocke, MM-SHR] - S-q, S-Rocke, MM-SHR are set to max from error sim
rngFun    = @(~) GaussDist.rand(p, n, muFun(p), sigmaFun(p));
b         = 0.5-0.5*p/n; % Maximum breakdown point
q         = 0.909;
gam       = 0.870;
mmshrCnst = 0.778;
stability_results{end+1} = stability_sim(               p, n, ntrials, epsilon, k, rngFun, mleFun, b, sqFun, sqDistParms, q, gam, mmshrCnst, tol, initFrac);
comp_eff_results{end+1}  = computational_efficiency_sim(p, n, ntrials, epsilon, k, rngFun,         b, sqFun, sqDistParms, q, gam, mmshrCnst, tol);


%% REPORT RESULTS
tf = now;

% Display results
for i = 1:length(stability_results)
    fprintf(1, '\n\nCOMPUTATIONAL SIMULATION RESULTS: p=%d, n=%d, epsilon=%.2f\n', stability_results{i}.test_details.p, stability_results{i}.test_details.n, stability_results{i}.test_details.epsilon);
    fprintf(1, '%12s%18s%18s\n', 'Estimator', 'Mean Divergence', 'Median Iter./run');
    estimators1 = fieldnames(stability_results{i}); estimators2 = fieldnames(comp_eff_results{i});
    if ~isequal(estimators1,estimators2); error('Unexpected Results'); end
    for k = 1:length(estimators1)
        if strcmp(estimators1{k}, 'test_details')
            continue;
        end
        fprintf(1, '%12s%18.1E%18.1E\n', estimators1{k}, stability_results{i}.(estimators1{k}).mean_SP_dkl_fhat_fHatMLE, comp_eff_results{i}.(estimators1{k}).median_iterations2convergence);
    end
end
fprintf(1, 'Total Time = %.2f minutes\n\n', (tf-t0)*24*60);

% Print Results
fid = fopen(result_file, 'w');
for i = 1:length(stability_results)
    fprintf(fid, 'COMPUTATIONAL SIMULATION RESULTS: p=%d, n=%d, epsilon=%.2f\n', stability_results{i}.test_details.p, stability_results{i}.test_details.n, stability_results{i}.test_details.epsilon);
    fprintf(fid, '%12s%18s%18s\n', 'Estimator', 'Mean Divergence', 'Median Iter./run');
    estimators1 = fieldnames(stability_results{i}); estimators2 = fieldnames(comp_eff_results{i});
    if ~isequal(estimators1,estimators2); error('Unexpected Results'); end
    for k = 1:length(estimators1)
        if strcmp(estimators1{k}, 'test_details')
            continue;
        end
        fprintf(fid, '%12s%18.1E%18.1E\n', estimators1{k}, stability_results{i}.(estimators1{k}).mean_SP_dkl_fhat_fHatMLE, comp_eff_results{i}.(estimators1{k}).median_iterations2convergence);
    end
end
fprintf(fid, 'Total Time = %.2f minutes\n\n', (tf-t0)*24*60);
fclose(fid);
end
