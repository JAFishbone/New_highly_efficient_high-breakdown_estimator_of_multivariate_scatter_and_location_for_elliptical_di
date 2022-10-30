% Copyright (C) 2020-2022 Justin A Fishbone
function test_results = stability_sim(p, n, ntrials, epsilon, k, rngFun, mleFun, b, sqFun, sqDistParms, q, gam, mmshrCnst, tol, initFrac)
% test_results = stability_sim(p, n, ntrials, epsilon, k, rngFun, mleFun, b, sqFun, sqDistParms, q, gam, mmshrCnst, tol, initFrac)
% Execute the estimator stability simulation
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

% Do not save the location and scatter results of every iteration
SAVE_INDV_RUN_RESULTS = false;

% Run simulation
tic
test_results0 = execute_sim(p, n, ntrials, epsilon, k, rngFun, mleFun, b, sqFun, sqDistParms, q, gam, mmshrCnst, tol, initFrac);
test_results = calc_results(test_results0, ntrials, p, SAVE_INDV_RUN_RESULTS);
test_results = save_test_metadata(test_results, p, n, ntrials, epsilon, k, rngFun, mleFun, b, sqFun, sqDistParms, q, gam, mmshrCnst, tol, initFrac, SAVE_INDV_RUN_RESULTS);

% Save results
save(sprintf('stability_sim-p%d-n%d-eps%d-initFrac%d.mat',p,n,round(100*epsilon),round(100*initFrac)), 'test_results','p','n','epsilon','initFrac');

% Display results
fprintf(1, '\n\nSTABILITY RESULTS: p=%d, n=%d, epsilon=%.2f\n', p, n, epsilon);
toc
fprintf(1, '%12s%18s\n', 'Estimator', 'Mean Divergence');
estimators = fieldnames(test_results);
for k = 1:length(estimators)
    if strcmp(estimators{k}, 'test_details')
        continue;
    end
    fprintf(1, '%12s%18.1E\n', estimators{k}, test_results.(estimators{k}).mean_SP_dkl_fhat_fHatMLE);
end

end

function test_results = execute_sim(p, n, ntrials, epsilon, k, rngFun, mleFun, b, sqFun, sqDistParms, q, gam, mmshrCnst, tol, initFrac)

% Initialize result variables
SRocke.muHat = zeros(p,1,ntrials,2); SRocke.sigmaHat = zeros(p,p,ntrials,2); 
mmSHR.muHat  = zeros(p,1,ntrials,2); mmSHR.sigmaHat  = zeros(p,p,ntrials,2); 
Sq.muHat     = zeros(p,1,ntrials,2); Sq.sigmaHat     = zeros(p,p,ntrials,2);

n_outliers = floor(epsilon*n); % n = nX * p
n_init     = floor(initFrac*n);

h = waitbar(0, 'Wait');
t0 = toc;
% Run tests
for i = 1:ntrials

    % Generate X -- Use cell array for different outlier values for each estimator
    X_clean = rngFun();
    X = repmat({X_clean}, 1, 3);
    if epsilon > 0
        for j = 1:3
            X{j}(1,1:n_outliers) = k(j); % [S-q, S-Rocke, MM-SHR]
        end
    end
    
    % Various initial estimates
    [mu0MLE,     sigma0MLE]     = mleFun(X_clean);
    [mu0FracMLE, sigma0FracMLE] = mleFun(X_clean(:, end-n_init+1:end));
    mu0    = {   mu0MLE,    mu0FracMLE};
    sigma0 = {sigma0MLE, sigma0FracMLE};
    
    % Run estimator for each initial estimate: j=1 -- MLE with n uncontaminated samples, j=2 -- MLE with initFrac*n uncontaminated samples
    for j = 1:2
        % S-q Estimator
        try       [Sq.muHat(:,1,i,j),           Sq.sigmaHat(:,:,i,j)] = sEst( X{1}, sqFun, sqDistParms,   q, mu0{j}, sigma0{j}, tol, b);
        catch;     Sq.muHat(:,1,i,j) = NaN;     Sq.sigmaHat(:,:,i,j) = NaN; end
        % S-Rocke
        try   [SRocke.muHat(:,1,i,j),       SRocke.sigmaHat(:,:,i,j)] = sEst( X{2}, 'ROCKE',        [], gam, mu0{j}, sigma0{j}, tol, b);
        catch; SRocke.muHat(:,1,i,j) = NaN; SRocke.sigmaHat(:,:,i,j) = NaN; end
        % SHR-MM
        try    [mmSHR.muHat(:,1,i,j),        mmSHR.sigmaHat(:,:,i,j)] = mmshr(X{3},               mmshrCnst, mu0{j}, sigma0{j}, tol, b);
        catch;  mmSHR.muHat(:,1,i,j) = NaN;  mmSHR.sigmaHat(:,:,i,j) = NaN; end
    end
    
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


function test_results = calc_results(test_results, ntrials, p, save_indv_run_results)
% Generate metrics for each estimator
estimators = fieldnames(test_results);
for k = 1:length(estimators)
    % Generate metrics for each Monte Carlo run
    for i = 1:ntrials
        % Make sure we have shape matrices
        sp_hat_mle = test_results.(estimators{k}).sigmaHat(:,:,i,1) / det(test_results.(estimators{k}).sigmaHat(:,:,i,1))^(1/p);
        sp_hat = test_results.(estimators{k}).sigmaHat(:,:,i,2) / det(test_results.(estimators{k}).sigmaHat(:,:,i,2))^(1/p);
        % D(fHatMLE, fHat)
        test_results.(estimators{k}).SP_dkl_fhat_fHatMLE(i) = ( trace(sp_hat_mle\sp_hat) - p );
    end
    
    % Generate summary metrics
    test_results.(estimators{k}).mean_SP_dkl_fhat_fHatMLE = nanmean(test_results.(estimators{k}).SP_dkl_fhat_fHatMLE);
    
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
