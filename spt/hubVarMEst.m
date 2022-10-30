% Copyright (C) 2020-2022 Justin A Fishbone

function sig2 = hubVarMEst(x)
    p = 0.98;
    c = chi2inv(p,1);
    b = integral(@(t) chi2pdf(t,1).*weight(t,c).*t, -inf, inf); % Gaussian consistency

    mu = median(x);
    r2 = (x - mu).^2;
    sig2 = mad(x)^2;
    sig2Prev = inf;
    err = inf;
    while err > 1e-12
        wSig = weight(r2/sig2,c);
        sig2 = mean(wSig.*r2)/b;
        err = abs(sig2 - sig2Prev)/sig2;
        sig2Prev = sig2;
    end
end

function w = weight(t,c)
    w = min(1, c./t);
end