% Copyright (C) 2020-2022 Justin A Fishbone
function efficiencies = calc_efficiencies(p, beta, fd, dphi, sqFun, sqFixed, wMLE, maxInt, q, gam, mmPrm)

% S-q
tun = q;
sqFun.pEff = p;
sqFun.tun = tun;
rhoFun = @(t) sqFun.rhoFun(t);
wFun = @(t) sqFun.wFun(t);
sigma = fzero( @(sigma) integral( @(d) fd(d).*rhoFun(d/sigma), 0, inf ) - 0.5, [1e-6, 1e3]);
lambda1 = -2*beta* integral( @(d) 1/sigma * d.^(p/2+1) .* wFun(d/sigma) .* dphi(d), 0, inf );
zetaSq = lambda1^-2*p*(p+2) * integral( @(d) fd(d) .* (d/sigma).^2 .* wFun(d/sigma).^2, 0, inf );


% S-q MAX
tun = 1;
sqFun.pEff = p;
sqFun.tun = tun;
% rhoFun = @(t) sqFun.rhoFun(t);
wFun = @(t) sqFun.wFun(t);
sigma = 1;%fzero( @(sigma) integral( @(d) fd(d).*rhoFun(d/sigma), 0, inf ) - 0.5, [1e-6, 1e3]);
lambda1 = -2*beta* integral( @(d) 1/sigma * d.^(p/2+1) .* wFun(d/sigma) .* dphi(d), 0, inf );
zetaSqMax = lambda1^-2*p*(p+2) * integral( @(d) fd(d) .* (d/sigma).^2 .* wFun(d/sigma).^2, 0, inf );


% S-q FIXED
tun = q;
sqFixed.pEff = p;
sqFixed.tun = tun;
rhoFun = @(t) sqFixed.rhoFun(t);
wFun = @(t) sqFixed.wFun(t);
sigma = fzero( @(sigma) integral( @(d) fd(d).*rhoFun(d/sigma), 0, inf ) - 0.5, [1e-6, 1e5]);
lambda1 = -2*beta* integral( @(d) 1/sigma * d.^(p/2+1) .* wFun(d/sigma) .* dphi(d), 0, inf );
zetaSqFixed = lambda1^-2*p*(p+2) * integral( @(d) fd(d) .* (d/sigma).^2 .* wFun(d/sigma).^2, 0, inf );


% S-Rocke
tun = gam;
sFun = sRocke;
sFun.pEff = p;
sFun.tun = tun;
rhoFun = @(t) sFun.rhoFun(t);
wFun = @(t) sFun.wFun(t);
sigma = fzero( @(sigma) integral( @(d) fd(d).*rhoFun(d/sigma), 0, inf ) - 0.5, [1e-6, 1e8]);
lambda1 = -2*beta* integral( @(d) 1/sigma * d.^(p/2+1) .* wFun(d/sigma) .* dphi(d), 0, inf );
zetaRocke = lambda1^-2*p*(p+2) * integral( @(d) fd(d) .* (d/sigma).^2 .* wFun(d/sigma).^2, 0, inf );


% S-Bisq
sFun = sBisq;
sFun.pEff = p;
sFun.tun = tun;
rhoFun = @(t) sFun.rhoFun(t);
wFun = @(t) sFun.wFun(t);
sigma = fzero( @(sigma) integral( @(d) fd(d).*rhoFun(d/sigma), 0, inf ) - 0.5, [1e-6, 1e8]);
lambda1 = -2*beta* integral( @(d) 1/sigma * d.^(p/2+1) .* wFun(d/sigma) .* dphi(d), 0, inf );
zetaBisq = lambda1^-2*p*(p+2) * integral( @(d) fd(d) .* (d/sigma).^2 .* wFun(d/sigma).^2, 0, inf );


% MM-SHR
wMM = @(t) ( (-1.944+1.728*t-0.312*t.^2+0.016*t.^3) .* (t<=9 & t>4) + (t<=4) );
if isempty(mmPrm) % Find maximum efficiency
    lambda1 = @(sigma) -2*beta* integral( @(d) 1/sigma * d.^(p/2+1) .* wMM(d/sigma) .* dphi(d), 0, inf );
    zetaMMFun =  @(sigma) lambda1(sigma).^-2*p*(p+2) * integral( @(d) fd(d) .* (d/sigma).^2 .* wMM(d/sigma).^2, 0, inf );
    [tunSigmaa, zetaMMa] = fminbnd(zetaMMFun, 0.2, 1); %#ok<ASGLU>
    [tunSigmab, zetaMMb] = fminbnd(zetaMMFun, 1, 500); %#ok<ASGLU>
    [tunSigmac, zetaMMc] = fminbnd(zetaMMFun, 500, 1e4); %#ok<ASGLU>
    [tunSigmad, zetaMMd] = fminbnd(zetaMMFun, 1e4, 1e5); %#ok<ASGLU>
    [tunSigmae, zetaMMe] = fminbnd(zetaMMFun, 1e5, 1e7); %#ok<ASGLU>
    zetaMM = min([zetaMMa,zetaMMb,zetaMMc,zetaMMd,zetaMMe]);
else
    wMM = @(t) ( (-1.944+1.728*t-0.312*t.^2+0.016*t.^3) .* (t<=9 & t>4) + (t<=4) );
    dwMM = @(t) ( (1.728-0.624*t+0.048*t.^2) .* (t<=9 & t>4) );
    mmSigma = fzero( @(sigma) integral( @(d) fd(d).*mmRhoFun(d/sigma), 0, inf ) - 0.5, [1e-3, 1e8]);
    sigma = mmPrm*mmSigma;
    lambda1  = integral( @(d) fd(d) .* ((p+2)*(d/sigma).*wMM(d/sigma) + 2*(d/sigma).^2.*dwMM(d/sigma)), 0, inf);
    zetaMM = lambda1^-2*p*(p+2) * integral( @(d) fd(d) .* (d/sigma).^2 .* wMM(d/sigma).^2, 0, inf );
end


% SCM
zetaSample = p/(p+2) * integral( @(d) fd(d).*d.^2, 0, maxInt) /  integral( @(d) fd(d).*d, 0, maxInt).^2;

% M-MLE
psi = @(d) d.*wMLE(d);
DEL = 1e-7;
dFun = @(x, fun) (fun(x+DEL) - fun(x-DEL)) ./ (2*DEL);
dpsi = @(d) dFun(d, psi);
sigma = 1;
psi1 = integral( @(d) fd(d).*psi(d/sigma).^2, 2*DEL, maxInt)/(p^2 + 2*p);
psi2 = integral( @(d) fd(d).*dpsi(d/sigma).*d, 2*DEL, maxInt )/(p*sigma);
zetaMLE = (p+2)^2*psi1/(p+2*psi2)^2;


% Efficiencies
efficiencies.sq      = zetaMLE / zetaSq;
efficiencies.sqMax   = zetaMLE / zetaSqMax;
efficiencies.sqFixed = zetaMLE / zetaSqFixed;
efficiencies.sRocke  = zetaMLE / zetaRocke;
efficiencies.mmShr   = zetaMLE / zetaMM;
efficiencies.sBisq   = zetaMLE / zetaBisq;
efficiencies.sample  = zetaMLE / zetaSample;
end