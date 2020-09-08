function D2 = mcd(ZI, ZJ)
% mcd - mel-cepstral distortion
%
% Input:
%   ZI [1, n]
%   ZJ [m2, n]
% Output:
%   D2 [m2, 1]: distance between observations ZI and ZJ

const = 10 * sqrt(2) / log(10);
D2 = const * vecnorm(ZI - ZJ, 2, 2);
