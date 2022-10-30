% Copyright (C) 2020-2022 Justin A Fishbone
function test_results = calc_sim_results(test_results, ntrials, nsteps, ~, sigma, p, save_indv_run_results)
% test_results = calc_sim_results(test_results, ntrials, nsteps, mu, sigma, p, save_indv_run_results)
% Utility to compile simulation results for simulations
%

% Get estimators, and ensure MLE is 1st
estimators = fieldnames(test_results);
estimators(strcmp(estimators, 'ml')) = [];
estimators = ['ml'; estimators];

% For xi variance estimation
omega = sigma/det(sigma)^(1/p);
denom = triu((eye(p^2) + Kmn(p,p))*kron(omega,omega) - 2/p*omega(:)*omega(:)');

% Generate metrics for each estimator
for k = 1:length(estimators)
    for j = 1:nsteps
        % Generate metrics for each Monte Carlo run
        for i = 1:ntrials
            omegaHat = test_results.(estimators{k}).sigmaHat(:,:,i,j) / det(test_results.(estimators{k}).sigmaHat(:,:,i,j))^(1/p);
            
            % Metrics
            % D(fHat||f)
            test_results.(estimators{k}).SP_dkl_fhat_f(i,j) = ( trace(omega\omegaHat) - p );
            % H^2
            test_results.(estimators{k}).SP_h2(i,j) = 2 - det(omega*omegaHat)^(1/4)*2^((p+2)/2)/det(omega+omegaHat)^(1/2);
            % xi hat
            numRM(:,:,i) = (omegaHat(:) - omega(:))*(omegaHat(:) - omega(:))'; %#ok<AGROW>
        end
    % xi hat
    nn = sum( ~isnan( squeeze( test_results.(estimators{k}).sigmaHat(1,1,:,j) ) ) );
    numR = nn*nanmean(numRM,3);
    test_results.(estimators{k}).xiHatMatR(:,j) = numR(denom(:)~=0) ./ denom(denom(:)~=0);
    test_results.(estimators{k}).xiHatR(j)      = nanmean(test_results.(estimators{k}).xiHatMatR(:,j));
    end
    
    % Generate summary metrics
    % Mean   
    test_results.(estimators{k}).mean_SP_dkl_fhat_f = nanmean(test_results.(estimators{k}).SP_dkl_fhat_f, 1);
    test_results.(estimators{k}).mean_SP_h2         = nanmean(test_results.(estimators{k}).SP_h2, 1);
    % Efficiency
    test_results.(estimators{k}).eff_SP_dkl_fhat_f = test_results.ml.mean_SP_dkl_fhat_f ./ test_results.(estimators{k}).mean_SP_dkl_fhat_f;
    test_results.(estimators{k}).eff_SP_h2         = test_results.ml.mean_SP_h2 ./ test_results.(estimators{k}).mean_SP_h2;
    test_results.(estimators{k}).eff_SP_xi         = test_results.ml.xiHatR ./ test_results.(estimators{k}).xiHatR;

    % Remove individual run results to save memory/disk
    if ~save_indv_run_results
        test_results.(estimators{k}).muHat    = [];
        test_results.(estimators{k}).sigmaHat = [];
    end
end
end