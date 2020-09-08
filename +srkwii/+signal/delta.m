function y = delta(x, varargin)
% delta - delta calculation
%   y = delta(x, coeffs)
%   y = delta(x, order, width)

narginchk(2, 3);

parser = inputParser;
parser.addRequired('x');

if nargin == 2
  parser.addRequired('coeffs');
  parser.parse(x, varargin{:});
  c = varargin{1};
  c = c(:)';
elseif nargin == 3
  parser.addRequired('order', @(a) isscalar(a) && (a == 1 || a == 2));
  parser.addRequired('width', @(a) isscalar(a) && a >= 1);
  parser.parse(x, varargin{:});
  n = varargin{1};
  w = varargin{2};
  a1 = sum((1:w) .^ 2) * 2;
  if n == 1
    c = (-w:w) ./ a1;
  elseif n == 2
    a2 = sum((1:w) .^ 4) * 2;
    a0 = w * 2 + 1;
    c = 2 * (a0 * (-w:w) .^ 2 - a1) ./ (a2 * a0 - a1 * a1);
  end
end

l = numel(c);
T = size(x, 2);
x = [repmat(x(:, 1), 1, floor(l / 2)) x repmat(x(:, end), 1, floor(l / 2))];
y = cell2mat(arrayfun(@(n) sum(x(:, n:n + l - 1) .* c, 2), 1:T, 'un', 0));
