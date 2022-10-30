% Copyright (C) 2020-2022 Justin A Fishbone
function [muHat, omegaHat, err] = mmshr(X, const, mu0, omega0, tol, b)
% [muHat, omegaHat, err] = mmshr(X, const, mu0, omega0, tol, b)
% MM-Estimator with Smooth Hard Rejection (SHR) rho function
%

% Get dimension and sample size
p = size(X,1);
n = size(X,2);

% Determine b
if nargin < 6 || isempty(b)
    b = 0.5*(1-p/n);
end

% Convergence Tolerance
if nargin < 5 || isempty(tol)
    tol = 1e-6;
end

% Initial values for mu, omega, and squared Mahalanobis distances
muPrev = mu0;
omegaPrev = omega0/det(omega0)^(1/p);
d = real(diag((X-muPrev)'*(omegaPrev\(X-muPrev))).');

% Solve for S using Eq. 6.44 in Maronna et al. 2009
S0 = median(d)/6.5;
lb = S0;
while (mean(rhoFun(d/lb)) - b) < 0
    lb = lb/2;
end
ub = S0;
while (mean(rhoFun(d/ub)) - b) > 0 && ub ~= inf
    ub = ub*2;
end
solRng = [lb, ub];
S = fzero( (@(s)mean(rhoFun(d./s)) - b), solRng );

% Loop until converged
err(1) = inf; j = 1; % Will use length of err to keep track of how many iterations until convergence
tmp = zeros(p,p,n);
while err(j) > tol
    % Weights
    w = wFun(d/(const*S));

    % Estimate mu
    muNext = sum(w.*X,2)/sum(w);
    % Estimate sigma
    for i = 1:n
        tmp(:,:,i) = (X(:,i)-muNext)*(X(:,i)-muNext)';
    end
    tmp2 = reshape(w,1,1,n) .* tmp;
    omegaNext = mean(tmp2,3);
    while det(omegaNext) > 1e100 % Sometimes omegaNext can get too large and det(omegaNext) computes to inf in MATLAB, so shrink it
        omegaNext = omegaNext/10;
    end
    omegaNext = omegaNext/det(omegaNext)^(1/p);
    % For now, error metric is Gaussian KL-divergence between previous and next estimates
    if sum(isnan(omegaNext(:)))
        error('error');
    end
    j = j + 1;
    % err = trace(omegaNext\omegaPrev) - p - log(det(omegaPrev)/det(omegaNext)) + (muNext-muPrev)'*(omegaNext\(muNext-muPrev));
    err(j) = trace(omegaNext\omegaPrev) - p; %#ok<AGROW>

    % Prepare for next iteration
    muPrev    = muNext;
    omegaPrev = omegaNext;
    d = real(diag((X-muPrev)'*(omegaPrev\(X-muPrev))).');
end
muHat    = muNext;
omegaHat = omegaNext;
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% SUPPORT FUNCTIONS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function rho = rhoFun(t)
% Rho function is the antiderivative of the weight function

% Constants for normalizing rho to max=1
at4 = -1.944*4+1.728/2*4^2-0.312/3*4^3+0.016/4*4^4;
at9 = -1.944*9+1.728/2*9^2-0.312/3*9^3+0.016/4*9^4;
normCnst = at9 + 4 - at4;

rho = ( -1.944*t+1.728/2*t.^2-0.312/3*t.^3+0.016/4*t.^4 + (4-at4) ) / normCnst ;
rho(t<=4) = t(t<=4)/normCnst;
rho(t>9)  = 1;
end

function w = wFun(t)
% Equation 6.47 from Maronna et al. 2019
w       = -1.944+1.728*t-0.312*t.^2+0.016*t.^3;
w(t<=4) = 1;
w(t>9)  = 0;
end