function D2 = kl(ZI, ZJ)
% kl - Kullback-Leibler (KL) divergence
%
% Input:
%   ZI [1, n]
%   ZJ [m2, n]
% Output:
%   D2 [m2, 1]: distance between observations ZI and ZJ

D2 = sum(ZI .* (log(ZI) - log(ZJ)), 2);
