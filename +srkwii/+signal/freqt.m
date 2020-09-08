function y = freqt(input, varargin)
% freqt - frequency transformation
%   output = freqt(input)
%   output = freqt(___, Name, Value)
%
% Options:
%   OutputOrder: order of warped sequence (default: same as input)
%   InputAlpha : all-pass constant of input sequence
%   OutputAlpha: all-pass constant of output sequence
%
% References:
%   K. Tokuda, et al. Recursive Calculation of Mel-Cepstrum from LP Coefficients. Technical Report of Nagoya Institute of Technology, 1994.

parser = inputParser;
parser.addRequired('input');
parser.addOptional('OutputOrder', size(input, 1) - 1, @(a) isscalar(a) && 0 <= a);
parser.addOptional('InputAlpha', 0, @(a) isscalar(a) && -1 < a && a < 1);
parser.addOptional('OutputAlpha', 0.35, @(a) isscalar(a) && -1 < a && a < 1);
parser.parse(input, varargin{:});

a1 = parser.Results.InputAlpha;
a2 = parser.Results.OutputAlpha;
a = (a2 - a1) / (1 - a1 * a2);
b = 1 - a * a;

m2 = parser.Results.OutputOrder + 1;
x = input;
y = zeros(m2, size(x, 2));

for m = size(x, 1):-1:1
  z = y(2, :);
  y(2, :) = b * y(1, :) + a * y(2, :);
  y(1, :) = x(m, :) + a * y(1, :);
  for n = 3:m2
    t = z + a * (y(n, :) - y(n - 1, :));
    z = y(n, :);
    y(n, :) = t;
  end
end
