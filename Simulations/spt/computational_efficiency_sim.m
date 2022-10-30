% Copyright (C) 2020-2022 Justin A Fishbone
function test_results = computational_efficiency_sim(p, n, ntrials, epsilon, k, rngFun, b, sqFun, sqDistParms, q, gam, mmshrCnst, tol)
% test_results = computational_efficiency_sim(p, n, ntrials, epsilon, k, rngFun, b, sqFun, sqDistParms, q, gam, mmshrCnst, tol)
% Execute the estimator computational efficiency simulation
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

% Do not save the location and scatter results of every iteration
SAVE_INDV_RUN_RESULTS = false;

% Run simulation
tic
test_results0 = execute_sim(p, n, ntrials, epsilon, k, rngFun, b, sqFun, sqDistParms, q, gam, mmshrCnst, tol);
test_results = calc_results(test_results0, SAVE_INDV_RUN_RESULTS);
test_results = save_test_metadata(test_results, p, n, ntrials, epsilon, k, rngFun, b, sqFun, sqDistParms, q, gam, mmshrCnst, tol, SAVE_INDV_RUN_RESULTS);

% Save results
save(sprintf('comp_efficiency_sim-p%d-n%d-eps%d.mat',p,n,round(100*epsilon)), 'test_results','p','n','epsilon');

% Display results
fprintf(1, '\n\nCOMPUTATIONAL EFFICIENCY RESULTS: p=%d, n=%d, epsilon=%.2f\n', p, n, epsilon);
toc
fprintf(1, '%12s%18s\n', 'Estimator', 'Median Iter./Run');
estimators = fieldnames(test_results);
for k = 1:length(estimators)
    if strcmp(estimators{k}, 'test_details')
        continue;
    end
    fprintf(1, '%12s%18E\n', estimators{k}, test_results.(estimators{k}).median_iterations2convergence);
end
end


function test_results = execute_sim(p, n, ntrials, epsilon, k, rngFun, b, sqFun, sqDistParms, q, gam, mmshrCnst, tol)

% Initialize result variables
Sq.muHat     = zeros(p,1,ntrials); Sq.sigmaHat     = zeros(p,p,ntrials); Sq.num2convrg     = zeros(1,1,ntrials);
SRocke.muHat = zeros(p,1,ntrials); SRocke.sigmaHat = zeros(p,p,ntrials); SRocke.num2convrg = zeros(1,1,ntrials);
mmSHR.muHat  = zeros(p,1,ntrials); mmSHR.sigmaHat  = zeros(p,p,ntrials); mmSHR.num2convrg  = zeros(1,1,ntrials);

n_outliers = floor(epsilon*n);

% Run tests
h = waitbar(0, 'Wait');
t0 = toc;
for i = 1:ntrials
    % Generate X -- Use cell array for different outlier values for each estimator
    X_clean = rngFun();
    X = repmat({X_clean}, 1, 3);
    if epsilon > 0
        for j = 1:3
            X{j}(1,1:n_outliers) = k(j); % [S-q, S-Rocke, MM-SHR]
        end
    end
    
    % Determine initial estimates 
    [~, ~, m, sigma0] = KurtSDFull(X{1}.');
    mu0 = m.';
    % S-q Estimator
    try       [Sq.muHat(:,1,i),           Sq.sigmaHat(:,:,i), ~, sigHat] = sEst( X{1},   sqFun, sqDistParms,   q, mu0, sigma0, tol, b);     Sq.num2convrg(1,1,i) = length(sigHat);
    catch;     Sq.muHat(:,1,i) = NaN;     Sq.sigmaHat(:,:,i) = NaN;     Sq.num2convrg(1,1,i) = NaN; end

    % Determine initial estimates 
    [~, ~, m, sigma0] = KurtSDFull(X{2}.');
    mu0 = m.';
    % S-Rocke
    try   [SRocke.muHat(:,1,i),       SRocke.sigmaHat(:,:,i), ~, sigHat] = sEst( X{2}, 'ROCKE',          [], gam, mu0, sigma0, tol, b); SRocke.num2convrg(1,1,i) = length(sigHat);
    catch; SRocke.muHat(:,1,i) = NaN; SRocke.sigmaHat(:,:,i) = NaN; SRocke.num2convrg(1,1,i) = NaN; end
    
    
    % Determine initial estimates 
    [~, ~, m, sigma0] = KurtSDFull(X{3}.');
    mu0 = m.';
    % MM-SHR
    try    [mmSHR.muHat(:,1,i),        mmSHR.sigmaHat(:,:,i), err]       = mmshr(X{3},                 mmshrCnst, mu0, sigma0, tol, b);  mmSHR.num2convrg(1,1,i) = length(err);
    catch;  mmSHR.muHat(:,1,i) = NaN;  mmSHR.sigmaHat(:,:,i) = NaN;  mmSHR.num2convrg(1,1,i) = NaN; end
    
    % Update waitbar
    t = ntrials/i*(toc-t0)/60; tl = t-(toc-t0)/60;
    waitbar(i/ntrials, h, sprintf('Wait: %.1f min remaining of %.1f minutes total', tl, t));
end
close(h);

% Pack results
test_results.Sq     = Sq;
test_results.SRocke = SRocke;
test_results.mmSHR  = mmSHR;
end


function test_results = calc_results(test_results, save_indv_run_results)
% Get estimators, and ensure MLE is 1st

% Generate metrics for each estimator
estimators = fieldnames(test_results);
for k = 1:length(estimators)
    % Generate summary metrics
    test_results.(estimators{k}).median_iterations2convergence = nanmedian(test_results.(estimators{k}).num2convrg);

    % Remove individual run results to save memory/disk
    if ~save_indv_run_results
        test_results.(estimators{k}).muHat    = [];
        test_results.(estimators{k}).sigmaHat = [];
    end
end
end


function test_results = save_test_metadata(test_results, varargin)
for i = 2:nargin
    test_results.test_details.(inputname(i)) = varargin{i-1};
end
end