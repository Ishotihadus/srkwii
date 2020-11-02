function x = griffinlim(s, varargin)
% griffinlim - griffinlim phase reconstruction
%   x = griffinlim(s)
%   x = griffinlim(___, Name, Value)
%
% Options:
%   Window :
%   Period :
%   IterMax:
%   Alpha  :
%
% References:
%   D. W. Griffin and J. S. Lim. Signal estimation from modified short-time Fourier transform. IEEE Transactions on Acoustics, Speech, and Signal Processing, 32 (2), 1984.
%   N. Perraudin, et al. A fast Griffin-Lim algorithm. In Proc. 2013 IEEE Workshop on Applications of Signal Processing to Audio and Acoustics, 2013.

parser = inputParser;
parser.addRequired('input');
parser.addOptional('Window', blackman(256));
parser.addOptional('Period', 100, @(a) isscalar(a) && 1 <= a);
parser.addOptional('IterMax', 0, @(a) isscalar(a) && 0 <= a);
parser.addOptional('Alpha', 0, @(a) isscalar(a) && 0 <= a && a <= 1);
parser.parse(s, varargin{:});

window = parser.Results.Window(:);
period = parser.Results.Period;

M = size(s, 2);
L = (size(s, 1) - 1) * 2;

if isreal(s)
  phase = randn(L / 2 + 1, M);
else
  phase = angle(s);
end
s = abs(s);

% length is (period * (M - 1) + L)
den = conv([repmat([1; zeros(period - 1, 1)], M - 1, 1); 1], window .^ 2);

% avoid division by 0
den(den == 0) = 1;

if parser.Results.IterMax > 0 && parser.Results.Alpha ~= 0
  x = revert(s, phase, window, period, den);
  ybef = fft(cell2mat(arrayfun(@(m) x(m * period + 1:m * period + L), 0:M - 1, 'un', 0)) .* window);
  ybef = ybef(1:L / 2 + 1, :);
end

for iter = 1:parser.Results.IterMax
  x = revert(s, phase, window, period, den);
  y = fft(cell2mat(arrayfun(@(m) x(m * period + 1:m * period + L), 0:M - 1, 'un', 0)) .* window);
  y = y(1:L / 2 + 1, :);
  if parser.Results.Alpha ~= 0
    y = y + parser.Results.Alpha * (y - ybef);
    ybef = y;
  end
  phase = angle(y);
  fprintf('iteration #%02d - %g\n', iter, norm(abs(y) - s) ./ numel(y));
end

x = revert(s, phase, window, period, den);
x = x(floor(L / 2) + 1:end - ceil(L / 2));


function x = revert(s, phase, window, period, den)
s = s .* exp(1i * phase);
s = [s; conj(s(end-1:-1:2, :))];
s = real(ifft(s)) .* window;
x = zeros(size(den));
L = size(s, 1);
for m = 1:size(s, 2)
  pos = (m - 1) * period + (1:L);
  x(pos) = x(pos) + s(:, m);
end
x = x ./ den;
