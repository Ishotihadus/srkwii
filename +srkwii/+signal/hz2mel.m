function mel = hz2mel(hz)
% hz2mel - convert from hertz to mel scale
%   mel = hz2mel(hz)

mel = 2595 * log10(hz / 700 + 1);
