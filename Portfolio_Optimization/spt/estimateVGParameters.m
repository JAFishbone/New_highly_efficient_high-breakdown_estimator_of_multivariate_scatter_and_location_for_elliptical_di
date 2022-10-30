% Copyright (C) 2020-2022 Justin A Fishbone
function results = estimateVGParameters(X, q, mmTun, gam, mq)
% function results = estimateVGParameters(X, q, mmTun, alpha, mq)
% Estimate symmetric multivariate variance gamma distribution parameters
%

%%% Constants
% Convergence tolerance 
TOL = 1e-6;
b = []; % Use the default b (which is the maximum breakdown point) of the estimators

% Important values
% Dimension
p = size(X,1);

%%% Initial estimates using MLE, which appears convex based on testing of different ICs
mu0 = zeros(p,1);
omega0 = eye(p);
lambda0 = 1;
psi0 = 1;
%%% Pack initial estimates into structures for each estimator
ml.muHat{1} = mu0;
ml.omegaHat{1} = omega0;
ml.lambdaHat{1} = lambda0;
ml.psiHat{1}=psi0;
% Set initial guess in MLE's structure -- redundant variables, but simplifies code in estimateLambdaPsi() below
ml.prmsHat{1} = [lambda0, psi0];
%%% Iterative robust estimation of parameters for MLE and S-q
k = 2;
err = inf;
while err > TOL
    %%% Estimate location and shape
    [ml.muHat{k}, ~, ml.omegaHat{k}] = GenHyperDist.mle(X, ml.muHat{k-1}, ml.omegaHat{k-1}, TOL, ml.lambdaHat{k-1}, ml.psiHat{k-1}, 0);
    %%% Estimate lambda and psi
    d = diag(real((X-ml.muHat{k})'*(ml.omegaHat{k}\(X-ml.muHat{k}))));
    f = @(lmb,psi) psi^(p/4+lmb/2)/(2^(p/2+lmb-1)*gamma(lmb)*gamma(p/2)) * d.^(p/2-1) .* d.^(lmb/2-p/4) .* besselk(lmb-p/2,sqrt(psi*d));
    ml.prmsHat{k} = fmincon(@(prms) -sum(log(f(prms(1),prms(2)))), ml.prmsHat{k-1}, [],[],[],[],[0,0], [300,1e7],[],optimoptions('fmincon', 'Display', 'notify-detailed'));
    ml.lambdaHat{k} = ml.prmsHat{k}(1);
    ml.psiHat{k}    = ml.prmsHat{k}(2);
    % Check solution
    if ml.lambdaHat{k} > 300-1 || ml.lambdaHat{k} < 0.1
        error('ml.lambdaHat{%d} out of bounds',k);
    elseif ml.psiHat{k} > 1e7-1e6 || ml.psiHat{k} < 0.1
        error('ml.psiHat{%d} out of bounds',k);
    end
    %%% Calculate convergence metric as max of all metrics
    % For now, error metric of location and shape is Gaussian KL-divergence between last two iterations
    % Lambda and psi are proportion change
    err(k) = max([ abs(ml.psiHat{k}    - ml.psiHat{k-1})    /ml.psiHat{k}, ...
                   abs(ml.lambdaHat{k} - ml.lambdaHat{k-1}) /ml.lambdaHat{k}, ...
                   abs(trace(ml.omegaHat{k} \ml.omegaHat{k-1})  - p - (ml.muHat{k} -ml.muHat{k-1})' *(ml.omegaHat{k} \(ml.muHat{k}-ml.muHat{k-1}))) ]);
    if k == 50
        error('Max k');
    end
    k = k + 1;
end

%%% Use MLE as initial guess for robust estimators
mu0     = ml.muHat{end};
omega0  = ml.omegaHat{end};
lambda0 = ml.lambdaHat{end};
psi0    = ml.psiHat{end};

%%% Pack initial estimates into structures for each estimator
Sq.muHat{1}=mu0;
Sq.omegaHat{1}=omega0;
Sq.lambdaHat{1}=lambda0; SRocke.lambdaHat{1}=lambda0; mmSHR.lambdaHat{1}=lambda0;
Sq.psiHat{1}=psi0; SRocke.psiHat{1}=psi0; mmSHR.psiHat{1}=psi0; 

% Set initial guess in each estimator's structure -- redundant variables, but simplifies code in estimateLambdaPsi() below
Sq.prmsHat{1}=[lambda0, psi0]; SRocke.prmsHat{1}=[lambda0, psi0]; mmSHR.prmsHat{1}=[lambda0, psi0];
Sq.prmsHatBiased{1} = [lambda0, psi0]; SRocke.prmsHatBiased{1}=[lambda0, psi0]; mmSHR.prmsHatBiased{1}=[lambda0, psi0];

%%% Robust estimation of location and shape for S-Rocke, and MM-SHR estimators
% Rocke S-Estimator
[SRocke.muHat{1}, SRocke.omegaHat{1}] = sEst( X, 'ROCKE', [], gam, mu0, omega0, TOL, b);
SRocke = estimateLambdaPsi(X, p, mq, SRocke);
% MM-SHR Estimator
[mmSHR.muHat{1}, mmSHR.omegaHat{1}] = mmshr(X, mmTun, mu0, omega0, TOL, b);
mmSHR = estimateLambdaPsi(X, p, mq, mmSHR);

%%% Iterative robust estimation of parameters Sq
k = 2;
err = inf;
while err > TOL

    %%% Estimate location and shape
    sqDistParms.lambda = Sq.lambdaHat{k-1}; sqDistParms.psi = Sq.psiHat{k-1};
    [Sq.muHat{k}, Sq.omegaHat{k}] = sEst( X, 'QVARGAMMA', sqDistParms, q, Sq.muHat{k-1}, Sq.omegaHat{k-1}, TOL, b);
    Sq.omegaHat{k}   = Sq.omegaHat{k}/det(Sq.omegaHat{k})^(1/p);

    %%% Estimate lambda and psi
    Sq = estimateLambdaPsi(X, p, mq, Sq);        
    
    %%% Calculate convergence metric as max of all metrics
    % For now, error metric of location and shape is Gaussian KL-divergence between last two iterations
    % Lambda and psi are proportion change
    err(k) = max([ abs(Sq.psiHat{k}   - Sq.psiHat{k-1})   /Sq.psiHat{k}, ...
           0.001*  abs(Sq.lambdaHat{k}- Sq.lambdaHat{k-1})/Sq.lambdaHat{k}, ... % Adjust error criteria to account for steady-state changes
                   abs(trace(Sq.omegaHat{k}\Sq.omegaHat{k-1}) - p - (Sq.muHat{k}-Sq.muHat{k-1})'*(Sq.omegaHat{k}\(Sq.muHat{k}-Sq.muHat{k-1}))) ]);
    
    if k == 50
        error('Max k');
    end
    k = k + 1;
end

%%% Remove intermediate estimates and pack results, and add sample mean and covariance estimates
results.sample.muHat = mean(X,2);
results.sample.covMatHat = cov(X');
results.sample.omegaHat = results.sample.covMatHat/det(results.sample.covMatHat)^(1/p);
results.sq    = removeIntermediateResults(Sq);
results.SRocke = removeIntermediateResults(SRocke);
results.mmSHR = removeIntermediateResults(mmSHR);

%%% Scale omega to covariance matrix (assumes Sigma=Omega with scaling absorbed by psi)
results.sq.covMatHat    = covScale(results.sq.psiHat,    results.sq.lambdaHat)    * results.sq.omegaHat;
results.SRocke.covMatHat = covScale(results.SRocke.psiHat, results.SRocke.lambdaHat) * results.SRocke.omegaHat;
results.mmSHR.covMatHat = covScale(results.mmSHR.psiHat, results.mmSHR.lambdaHat) * results.mmSHR.omegaHat;

%%% Suppress warning messages about not pre-allocating arrays
%#ok<*AGROW>
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%% SUPPORT FUNCTIONS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function estStruct = estimateLambdaPsi(X, p, mq, estStruct)
% Estimate lambda and psi using M-Lq based on density of Mahalanobis distance.

d = diag(real((X-estStruct.muHat{end})'*(estStruct.omegaHat{end}\(X-estStruct.muHat{end}))));
f = @(lmb,psi) psi^(p/4+lmb/2)/(2^(p/2+lmb-1)*gamma(lmb)*gamma(p/2)) * d.^(p/2-1) .* d.^(lmb/2-p/4) .* besselk(lmb-p/2,sqrt(psi*d));
estStruct.prmsHatBiased{end+1} = fmincon(@(prms) -sum(f(prms(1),prms(2)).^(1-mq)), estStruct.prmsHatBiased{end}, [],[],[],[],[0,0], [100,1e7],[],optimoptions('fmincon', 'Display', 'notify-detailed'));
% Check solution
if estStruct.prmsHatBiased{end}(1) > 100-1
    error('lambdaHat out of bounds, %s, k=%d', inputname(4), length(estStruct.prmsHatBiased));
elseif estStruct.prmsHatBiased{end}(2) > 1e7-1e6 || estStruct.prmsHatBiased{end}(2) < 0.1
    error('psiHat out of bounds, %s, k=%d', inputname(4), length(estStruct.prmsHatBiased));
end
    
% Try to remove some asymptotic bias
f_d2 = @(d,prms) prms(2)^(p/4+prms(1)/2)/(2^(p/2+prms(1)-1)*gamma(prms(1))*gamma(p/2)) * d.^(p/2-1) .* d.^(prms(1)/2-p/4) .* besselk(prms(1)-p/2,sqrt(prms(2)*d));
estStruct.prmsHat{end+1} = ...
    fsolve( @(prmsHat) ...
        estStruct.prmsHatBiased{end} - ...
        fmincon( @(prmsBiasedAsymptotic) integral( @(d) -f_d2(d,prmsBiasedAsymptotic).^(1-mq) .* f_d2(d,prmsHat), 0, inf), estStruct.prmsHatBiased{end}, [],[],[],[],[0,0],[100,1e7],[],optimoptions('fmincon', 'Display', 'notify-detailed')),...
        estStruct.prmsHatBiased{end}, optimoptions('fsolve', 'Display', 'off'));
% Check solution
if 0.1 < sum(abs(estStruct.prmsHatBiased{end} - fmincon( @(prmsAsymptotic) integral( @(d) -f_d2(d,prmsAsymptotic).^(1-mq) .* f_d2(d,estStruct.prmsHat{end}), 0, inf), estStruct.prmsHatBiased{end}, [],[],[],[],[0,0],[100,1e7],[],optimoptions('fmincon', 'Display', 'notify-detailed'))))
    warning('Unable to remove all bias for %s, k=%d', inputname(4), length(estStruct.prmsHat));
end
estStruct.lambdaHat{end+1} = estStruct.prmsHat{end}(1);
estStruct.psiHat{end+1}    = estStruct.prmsHat{end}(2);
% Check solution
if estStruct.lambdaHat{end} > 100-1 || estStruct.lambdaHat{end} < 0.1
    error('lambdaHat out of bounds, %s, k=%d', inputname(4), length(estStruct.prmsHatBiased));
elseif estStruct.psiHat{end} > 1e7-1e6 || estStruct.psiHat{end} < 0.1
    error('psiHat out of bounds, %s, k=%d', inputname(4), length(estStruct.prmsHatBiased));
end

end


function outStruct = removeIntermediateResults(inStruct)
% Replace all cell arrays in resultStruct with the last element

fields = fieldnames(inStruct);
outStruct = inStruct;
for i = 1:length(fields)
    if iscell(outStruct.(fields{i}))
        outStruct.(fields{i}) = outStruct.(fields{i}){end};
    end
end
end


function scaleFactor = covScale(psi, lambda)
% Calculate scale factor to scale sigma to covariance matrix

% Since chi=0 and function is undefined at 0, take the limit by using a small chi
SM_CHI = 1e-10;

scaleFactor = sqrt(SM_CHI/psi)*besselk(lambda+1, sqrt(SM_CHI*psi))/besselk(lambda, sqrt(SM_CHI*psi));
end