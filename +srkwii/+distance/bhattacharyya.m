function D2 = bhattacharyya(ZI, ZJ)
% kl - Kullback-Leibler (KL) divergence
%
% Input:
%   ZI [1, n]
%   ZJ [m2, n]
% Output:
%   D2 [m2, 1]: distance between observations ZI and ZJ

D2 = -log(sum(sqrt(ZI .* ZJ), 2));
