% Copyright (C) 2020-2022 Justin A Fishbone
function K = Kmn(m,n)
% function K = Kmn(m,n)
% Return commutation matrix
K = zeros(m*n, m*n);
m0 = 1:(m*n);
N = reshape(m0, m,n)';
n0 = N(:);
for i = 1:(m*n)
   K(m0(i), n0(i)) = 1;
end
end
