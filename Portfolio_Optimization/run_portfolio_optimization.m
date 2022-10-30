% Copyright (C) 2020-2022 Justin A Fishbone
function run_portfolio_optimization
% Run financial portfolio optimizations


%% Configuration Parameters
DESIRED_DAILY_RETURNS = 1.1^(1/252)-1;
% Estimator tuning parameters
q     = 0.998; % S-q tuning parameter
mq    = 0.9999; % M-Lq tuning parameter
mmTun = 1.5; % MM-SHR tuning parameter
gam   = 1; % S-Rocke tuning parameter


%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% First run: 2020Q1 vs 2020Q1-Pre-pandemic
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

runPairs = { ...
'y20q1Inds', 'y20q1PreCrash'
};

% Run
for i = 1:size(runPairs,1)
    results1{i} = portfolio_optimization(runPairs{i,1}, runPairs{i,2}, DESIRED_DAILY_RETURNS, q, mq, mmTun, gam, true); %#ok<*AGROW>
end


%%% Print summary metrics
fprintf(1, '\n\nSUMMARY METRICS...\n');
fields = fieldnames(results1{1});
for i = 1:size(runPairs,1)
    fprintf(1, '\nPAIR: %s-%s\n', runPairs{i,1}, runPairs{i,2});
    % Variance of daily returns
    for j = 1:length(fields)
        fprintf(1, 'Variance of daily returns, %6s: %.1f\n', fields{j}, var(results1{i}.(fields{j}).dailyReturns)*1e6);
    end
    % Variance of daily returns
    for j = 1:length(fields)
        fprintf(1, 'Robust variance of daily returns, %6s: %.1f\n', fields{j}, hubVarMEst(results1{i}.(fields{j}).dailyReturns)*1e6);
    end
end


%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% Run 4 Years of Annual Same-Year Performance
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

runPairs = { ...
'y16Inds', 'y16Inds'
'y17Inds', 'y17Inds'
'y18Inds', 'y18Inds'
'y19Inds', 'y19Inds'
};

% Run
for i = 1:size(runPairs,1)
    results3{i} = portfolio_optimization(runPairs{i,1}, runPairs{i,2}, DESIRED_DAILY_RETURNS, q, mq, mmTun, gam, true); %#ok<*SAGROW>
end


%%% Print summary metrics
fprintf(1, '\n\nSUMMARY METRICS...\n');
fields = fieldnames(results3{1});
variances = nan(length(fields), size(runPairs,1));
rvariances = nan(length(fields), size(runPairs,1));
for i = 1:size(runPairs,1)
    fprintf(1, '\nPAIR: %s-%s\n', runPairs{i,1}, runPairs{i,2});
    % Normalized variance of daily returns -- Ignore sample estimator
    for j = 1:length(fields)
        if strcmpi(fields{j}, 'sample')
            continue;
        end
        variances(j,i) = var(results3{i}.(fields{j}).dailyReturns)*1e6;
        rvariances(j,i) = hubVarMEst(results3{i}.(fields{j}).dailyReturns)*1e6;
    end
    for j = 1:length(fields)
        if strcmpi(fields{j}, 'sample')
            continue;
        end
        fprintf(1, 'Variance of daily returns, %6s: %.1f\n', fields{j}, variances(j,i));
    end
    for j = 1:length(fields)
        if strcmpi(fields{j}, 'sample')
            continue;
        end
        fprintf(1, 'Robust variance of daily returns, %6s: %.1f\n', fields{j}, rvariances(j,i));
    end
end
fprintf(1, '\nAVERAGES:\n');
for j = 1:length(fields)
    if strcmpi(fields{j}, 'sample')
        continue;
    end
    fprintf(1, 'Average of variances of daily returns, %6s: %.1f\n', fields{j}, nanmean(variances(j,:)));
end
fprintf(1, '\nAVERAGES:\n');
for j = 1:length(fields)
    if strcmpi(fields{j}, 'sample')
        continue;
    end
    fprintf(1, 'Average of robust variances of daily returns, %6s: %.1f\n', fields{j}, nanmean(rvariances(j,:)));
end

end