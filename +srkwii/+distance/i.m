function D2 = i(ZI, ZJ)
% i - I divergence a.k.a. generalized Kullback-Leibler divergence
%
% Input:
%   ZI [1, n]
%   ZJ [m2, n]
% Output:
%   D2 [m2, 1]: distance between observations ZI and ZJ

D2 = sum(ZI .* (log(ZI) - log(ZJ) - 1) + ZJ, 2);
