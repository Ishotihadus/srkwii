function output = frame(input, varargin)
% frame - extract frame from data sequence
%   output = frame(input)
%   output = frame(___, Name, Value)
%
% Options:
%   Length (=L) : frame length (default: 256)
%   Period      : frame period (default: 100)
%   Centering   : centering input(1) (default: true)
% Output:
%   output [L, T]: framed series

parser = inputParser;
parser.addRequired('input', @isreal);
parser.addOptional('Length', 256, @(a) isscalar(a) && isreal(a) && 1 <= a);
parser.addOptional('Period', 100, @(a) isscalar(a) && isreal(a) && 1 <= a);
parser.addOptional('Centering', true, @(a) isscalar(a));
parser.parse(input, varargin{:});

L = parser.Results.Length;
P = parser.Results.Period;
len = numel(input);
if parser.Results.Centering
  input = [zeros(floor(L / 2), 1); input(:); zeros(ceil(L / 2), 1)];
else
  input = [input(:); zeros(L, 1)];
end
output = input((1:L)' + (0:floor((len - 1) / P)) * P);
