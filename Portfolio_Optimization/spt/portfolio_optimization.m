% Copyright (C) 2020-2022 Justin A Fishbone
function results = portfolio_optimization(estInds, perfInds, DESIRED_DAILY_RETURNS, q, mq, mmTun, gam, SUPRESS_OUTPUT)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% NOTE, parameter ML results from FY19:
% psiHat    = 4.5937e+04
% lambdaHat = 2.3581
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%% CONSTANTS
DATA_FILE = 'linear_returns.csv';
if nargin < 3
    DESIRED_DAILY_RETURNS = 1.1^(1/252)-1;
    % Estimator tuning parameters
    q     = 0.998; % S-q tuning parameter
    mq    = 0.9999; % M-Lq tuning parameter
    mmTun = 1.5; % MM-SHR tuning parameter
    gam   = 1; % S-Rocke tuning parameter
    SUPRESS_OUTPUT = false;
end

%%% Load data
fid = fopen(DATA_FILE);
firstline = textscan(fgetl(fid), '%s', 'Delimiter', ',');
stockTkrs = firstline{1}(2:end);
data = textscan(fid, '%f', 'Delimiter', ',');
data = reshape(data{1}, 1+length(stockTkrs), length(data{1})/(1+length(stockTkrs)));
dateNums = data(1,:);
linReturnsAll = data(2:end,:) - 1;
fclose(fid);


%%% Determine indices of interest
y16Inds       = datenum('01-JAN-2016') <= dateNums & dateNums <= datenum('31-DEC-2016');
y16q1Inds     = datenum('01-JAN-2016') <= dateNums & dateNums <= datenum('31-MAR-2016');
y16q2Inds     = datenum('01-APR-2016') <= dateNums & dateNums <= datenum('30-JUN-2016');
y16q3Inds     = datenum('01-JUL-2016') <= dateNums & dateNums <= datenum('30-SEP-2016');
y16q4Inds     = datenum('01-OCT-2016') <= dateNums & dateNums <= datenum('31-DEC-2016');

y17Inds       = datenum('01-JAN-2017') <= dateNums & dateNums <= datenum('31-DEC-2017');
y17q1Inds     = datenum('01-JAN-2017') <= dateNums & dateNums <= datenum('31-MAR-2017');
y17q2Inds     = datenum('01-APR-2017') <= dateNums & dateNums <= datenum('30-JUN-2017');
y17q3Inds     = datenum('01-JUL-2017') <= dateNums & dateNums <= datenum('30-SEP-2017');
y17q4Inds     = datenum('01-OCT-2017') <= dateNums & dateNums <= datenum('31-DEC-2017');

y18Inds       = datenum('01-JAN-2018') <= dateNums & dateNums <= datenum('31-DEC-2018');
y18q1Inds     = datenum('01-JAN-2018') <= dateNums & dateNums <= datenum('31-MAR-2018');
y18q2Inds     = datenum('01-APR-2018') <= dateNums & dateNums <= datenum('30-JUN-2018');
y18q3Inds     = datenum('01-JUL-2018') <= dateNums & dateNums <= datenum('30-SEP-2018');
y18q4Inds     = datenum('01-OCT-2018') <= dateNums & dateNums <= datenum('31-DEC-2018');

y19Inds       = datenum('01-JAN-2019') <= dateNums & dateNums <= datenum('31-DEC-2019');
y19q1Inds     = datenum('01-JAN-2019') <= dateNums & dateNums <= datenum('31-MAR-2019');
y19q2Inds     = datenum('01-APR-2019') <= dateNums & dateNums <= datenum('30-JUN-2019');
y19q3Inds     = datenum('01-JUL-2019') <= dateNums & dateNums <= datenum('30-SEP-2019');
y19q4Inds     = datenum('01-OCT-2019') <= dateNums & dateNums <= datenum('31-DEC-2019');


y20q1Inds     = datenum('01-JAN-2020') <= dateNums & dateNums <= datenum('31-MAR-2020');
y20q1PreCrash = datenum('01-JAN-2020') <= dateNums & dateNums <= datenum('20-FEB-2020');
y20q2Inds     = datenum('01-APR-2020') <= dateNums & dateNums <= datenum('30-JUN-2020');
y20q3Inds     = datenum('01-JUL-2020') <= dateNums & dateNums <= datenum('30-SEP-2020');

y16_19Inds    = datenum('01-JAN-2016') <= dateNums & dateNums <= datenum('31-DEC-2019');
y19q12Inds    = datenum('01-JAN-2019') <= dateNums & dateNums <= datenum('30-JUN-2019');
y19q34Inds    = datenum('01-JUL-2019') <= dateNums & dateNums <= datenum('31-DEC-2019');
y19q45Inds    = datenum('01-OCT-2019') <= dateNums & dateNums <= datenum('31-MAR-2020');
y20q23Inds    = datenum('01-APR-2020') <= dateNums & dateNums <= datenum('30-SEP-2020');


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% SPECIFY WHICH DATA TO USE
if nargin < 2
    % Data for estimation
    EST_INDS  = y20q1Inds;
    % Data for characterizing performance
    PERF_INDS = y20q1PreCrash;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
else
    EST_INDS  = eval(estInds);
    PERF_INDS = eval(perfInds);
end
data4est  = linReturnsAll(:,EST_INDS);
data4perf = linReturnsAll(:,PERF_INDS);


%%% Estimate variance gamma parameters with various estimators
results = estimateVGParameters(data4est, q, mmTun, gam, mq);
% Each field of the results structure will correspond to an estimator
fields = fieldnames(results);

%%%
%%% Calculate allocations
for i = 1:length(fields)
    results.(fields{i}).allocations = optimalAllocation(DESIRED_DAILY_RETURNS, results.(fields{i}).muHat, results.(fields{i}).covMatHat);
end

%%% Determine actual returns
for i = 1:length(fields)
    results.(fields{i}).dailyReturns = results.(fields{i}).allocations'*data4perf; % Assumes daily re-balance
end

if ~SUPRESS_OUTPUT
    %%% Summary metrics
    fprintf(1, 'Lin desired daily return: %.6f\n',    (DESIRED_DAILY_RETURNS));
    % variance of daily returns
    for i = 1:length(fields)
        fprintf(1, 'Variance of daily returns, %6s: %.6f\n', fields{i}, var(results.(fields{i}).dailyReturns));
    end
end

end % run_sim


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%% SUPPORT FUNCTIONS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function allocations = optimalAllocation(desiredLinReturn, mu, covMat)
one = ones(size(mu));
allocations = (mu'*(covMat\mu))*(covMat\one) - (one'*(covMat\mu))*(covMat\mu) + ...
    desiredLinReturn*( (one'*(covMat\one))*(covMat\mu) - (mu'*(covMat\one))*(covMat\one));
allocations = allocations/sum(allocations);
end