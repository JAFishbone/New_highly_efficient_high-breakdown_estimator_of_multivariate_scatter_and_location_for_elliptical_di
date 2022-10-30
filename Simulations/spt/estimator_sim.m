% Copyright (C) 2020-2022 Justin A Fishbone
function test_results = estimator_sim(p, n, nsteps, ntrials, distParm, mu, sigma, epsilon, k, rngFun, mleFun, b, sqFun, sqDistParms, q, gam, mmshrCnst, inclBisq, initializeTrue, tol, sq2Fun)
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
% sq2Fun         - Second Sq function
% 

% Do not save the location and scatter results of every iteration
SAVE_INDV_RUN_RESULTS = false;

% Ensure vectors are correct length, and vectorize scalar inputs
distParm    = proc_vector_inputs(distParm, nsteps);
k           = proc_vector_inputs(k, nsteps);
sqDistParms = proc_vector_inputs(sqDistParms, nsteps);
q           = proc_vector_inputs(q, nsteps);
gam         = proc_vector_inputs(gam, nsteps);
mmshrCnst   = proc_vector_inputs(mmshrCnst, nsteps);

% Run test
tic;
if nargin < 21
    test_results0 = execute_sim(p, n, nsteps, ntrials, distParm, mu, sigma, epsilon, k, rngFun, mleFun, b, sqFun, sqDistParms, q, gam, mmshrCnst, inclBisq, initializeTrue, tol);
else
    test_results0 = execute_sim(p, n, nsteps, ntrials, distParm, mu, sigma, epsilon, k, rngFun, mleFun, b, sqFun, sqDistParms, q, gam, mmshrCnst, inclBisq, initializeTrue, tol, sq2Fun);
end
test_results  = calc_sim_results(test_results0, ntrials, length(distParm), mu, sigma, p, SAVE_INDV_RUN_RESULTS);
test_results  = save_test_metadata(test_results, p, n, nsteps, ntrials, distParm, mu, sigma, epsilon, k, rngFun, mleFun, b, sqFun, sqDistParms, q, gam, mmshrCnst, inclBisq, initializeTrue, tol, SAVE_INDV_RUN_RESULTS);
toc
end

function test_results = execute_sim(p, n, nsteps, ntrials, distParm, mu, sigma, epsilon, k, rngFun, mleFun, b, sqFun, sqDistParms, q, gam, mmshrCnst, inclBisq, initializeTrue, tol, sq2Fun)

% Initialize result variables
ml.muHat     = zeros(p,1,ntrials,nsteps); ml.sigmaHat     = zeros(p,p,ntrials,nsteps); 
SRocke.muHat = zeros(p,1,ntrials,nsteps); SRocke.sigmaHat = zeros(p,p,ntrials,nsteps); 
mmSHR.muHat  = zeros(p,1,ntrials,nsteps); mmSHR.sigmaHat  = zeros(p,p,ntrials,nsteps); 
Sq.muHat     = zeros(p,1,ntrials,nsteps); Sq.sigmaHat     = zeros(p,p,ntrials,nsteps);
Sq2.muHat    = zeros(p,1,ntrials,nsteps); Sq2.sigmaHat    = zeros(p,p,ntrials,nsteps);
SBisq.muHat  = zeros(p,1,ntrials,nsteps); SBisq.sigmaHat  = zeros(p,p,ntrials,nsteps);

n_outliers = floor(epsilon*n);

% Run tests
h = waitbar(0, 'Wait');
t0 = toc;
for j = 1:nsteps        
    for i = 1:ntrials

        % Generate X
        X = rngFun(distParm(j));
        if epsilon > 0
            X(1,1:n_outliers) = k(j);
        end
        
        % Determine initial estimates
        if initializeTrue
            mu0 = mu; sigma0 = sigma;
        else
            [~, ~, m, sigma0] = KurtSDFull(X.');
            mu0 = m.';
        end
        
        % MLE
        try       [ml.muHat(:,1,i,j),           ml.sigmaHat(:,:,i,j)] = mleFun(X, distParm(j), mu0, sigma0);
        catch;     ml.muHat(:,1,i,j) = NaN;     ml.sigmaHat(:,:,i,j) = NaN; end
        % S-q Estimator
        try       [Sq.muHat(:,1,i,j),           Sq.sigmaHat(:,:,i,j)] = sEst( X, sqFun, sqDistParms(j),   q(j), mu0, sigma0, tol, b);
        catch;     Sq.muHat(:,1,i,j) = NaN;     Sq.sigmaHat(:,:,i,j) = NaN; end
        % Second S-q Estimator
        if nargin > 20
        try       [Sq2.muHat(:,1,i,j),          Sq2.sigmaHat(:,:,i,j)] = sEst( X, sq2Fun, {},   q(j), mu0, sigma0, tol, b);
        catch;     Sq2.muHat(:,1,i,j) = NaN;    Sq2.sigmaHat(:,:,i,j) = NaN; end
        end
        % S-Rocke
        try   [SRocke.muHat(:,1,i,j),       SRocke.sigmaHat(:,:,i,j)] = sEst( X, 'ROCKE',           [], gam(j), mu0, sigma0, tol, b);
        catch; SRocke.muHat(:,1,i,j) = NaN; SRocke.sigmaHat(:,:,i,j) = NaN; end
        % SHR-MM
        try    [mmSHR.muHat(:,1,i,j),        mmSHR.sigmaHat(:,:,i,j)] = mmshr(X,                  mmshrCnst(j), mu0, sigma0, tol, b);
        catch;  mmSHR.muHat(:,1,i,j) = NaN;  mmSHR.sigmaHat(:,:,i,j) = NaN; end
        % S-Bisq
        if inclBisq
        try    [SBisq.muHat(:,1,i,j),        SBisq.sigmaHat(:,:,i,j)] = sEst( X, 'BISQ',            [],     [], mu0, sigma0, tol, b);
        catch;  SBisq.muHat(:,1,i,j) = NaN;  SBisq.sigmaHat(:,:,i,j) = NaN; end
        end
    end
    % Update waitbar
    t = nsteps/j*(toc-t0)/60; tl = t-(toc-t0)/60;
    waitbar(j/nsteps, h, sprintf('Wait: %.1f min remaining of %.1f minutes total', tl, t));
end
close(h);

% Pack results
test_results.ml = ml;
test_results.Sq = Sq;
test_results.SRocke = SRocke;
test_results.mmSHR = mmSHR;
if inclBisq
    test_results.SBisq = SBisq;
end
if nargin > 20
    test_results.Sq2 = Sq2;
end
end

function param = proc_vector_inputs(param, nsteps)
% Verify length of param, and vectorize it if it is a scalar
switch length(param)
    case {0, nsteps}
        % Do nothing
    case 1
        param = repmat(param, [1, nsteps]);
    otherwise
        error('Invalid length for parameter %s', inputname(1));
end
end


function test_results = save_test_metadata(test_results, p, n, nsteps, ntrials, distParm, mu, sigma, epsilon, k, rngFun, mleFun, b, sqFun, sqDistParms, q, gam, mmshrCnst, inclBisq, initializeTrue, tol, save_indv_run_results)
test_results.test_details.p = p;
test_results.test_details.n = n;
test_results.test_details.nsteps = nsteps;
test_results.test_details.ntrials = ntrials;
test_results.test_details.distParm = distParm;
test_results.test_details.mu = mu;
test_results.test_details.sigma = sigma;
test_results.test_details.epsilon = epsilon;
test_results.test_details.k = k;
test_results.test_details.rngFun = rngFun;
test_results.test_details.mleFun = mleFun;
test_results.test_details.b = b;
test_results.test_details.sqFun = sqFun;
test_results.test_details.sqDistParms = sqDistParms;
test_results.test_details.q = q;
test_results.test_details.gam = gam;
test_results.test_details.mmshrCnst = mmshrCnst;
test_results.test_details.inclBisq = inclBisq;
test_results.test_details.initializeTrue = initializeTrue;
test_results.test_details.tol = tol;
test_results.test_details.save_indv_run_results = save_indv_run_results;
end