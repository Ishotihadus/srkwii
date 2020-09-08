function D2 = is(ZI, ZJ)
% is - Itakura-Saito divergence
%
% Input:
%   ZI [1, n]
%   ZJ [m2, n]
% Output:
%   D2 [m2, 1]: distance between observations ZI and ZJ

D2 = sum(ZI ./ ZJ - log(ZI) + log(ZJ) - 1, 2);
