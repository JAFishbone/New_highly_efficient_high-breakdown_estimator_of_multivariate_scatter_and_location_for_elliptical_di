% Copyright (C) 2020-2022 Justin A Fishbone
function rho = mmRhoFun(t)
% Rho function is the antiderivative of the weight function

% Constants for normalizing rho to max=1
at4 = -1.944*4+1.728/2*4^2-0.312/3*4^3+0.016/4*4^4;
at9 = -1.944*9+1.728/2*9^2-0.312/3*9^3+0.016/4*9^4;
normCnst = at9 + 4 - at4;

rho = ( -1.944*t+1.728/2*t.^2-0.312/3*t.^3+0.016/4*t.^4 + (4-at4) ) / normCnst ;
rho(t<=4) = t(t<=4)/normCnst;
rho(t>9)  = 1;
end