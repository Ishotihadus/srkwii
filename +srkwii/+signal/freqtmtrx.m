function y = freqtmtrx(varargin)
% freqtmtrx - frequency transformation matrix
%   result = freqtmtrx()
%   result = freqtmtrx(___, Name, Value)
%
% Options:
%   InputOrder  [M]: order of minimum phase sequence
%   OutputOrder [N]: order of warped sequence (default: same as InputOrder)
%   InputAlpha     : all-pass constant of input sequence
%   OutputAlpha    : all-pass constant of output sequence
% Output:
%   y [M+1, N+1]: transformation matrix
%
% References:
%   K. Tokuda, et al. Recursive Calculation of Mel-Cepstrum from LP Coefficients. Technical Report of Nagoya Institute of Technology, 1994.

parser = inputParser;
parser.addOptional('InputOrder', 25, @(a) isscalar(a) && 0 <= a);
parser.addOptional('OutputOrder', false, @(a) isscalar(a) && (0 <= a || a == false));
parser.addOptional('InputAlpha', 0, @(a) isscalar(a) && -1 < a && a < 1);
parser.addOptional('OutputAlpha', 0.35, @(a) isscalar(a) && -1 < a && a < 1);
parser.parse(varargin{:});

a1 = parser.Results.InputAlpha;
a2 = parser.Results.OutputAlpha;
m1 = parser.Results.InputOrder + 1;
if parser.Results.OutputOrder
  m2 = parser.Results.OutputOrder + 1;
else
  m2 = m1;
end

x = eye(m1);
y = zeros(m2, m1);
a = (a2 - a1) / (1 - a1 * a2);
b = 1 - a * a;

for m = m1:-1:1
  z = y(2, :);
  y(2, :) = b * y(1, :) + a * y(2, :);
  y(1, :) = x(m, :) + a * y(1, :);
  for n = 3:m2
    t = z + a * (y(n, :) - y(n - 1, :));
    z = y(n, :);
    y(n, :) = t;
  end
end
