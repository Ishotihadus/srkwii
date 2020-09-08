function y = freqtmtrx(varargin)
% freqtmtrx - frequency transformation matrix
%   result = freqtmtrx()
%   result = freqtmtrx(___, Name, Value)
%
% Options:
%   InputOrder: order of minimum phase sequence
%   OutputOrder: order of warped sequence (default: same as InputOrder)
%   InputAlpha : all-pass constant of input sequence
%   OutputAlpha: all-pass constant of output sequence

parser = inputParser;
parser.addOptional('InputOrder', 25, @(a) isscalar(a) && 0 <= a);
parser.addOptional('OutputOrder', false, @(a) isscalar(a) && (0 <= a || a == false));
parser.addOptional('InputAlpha', 0, @(a) isscalar(a) && -1 < a && a < 1);
parser.addOptional('OutputAlpha', 0.35, @(a) isscalar(a) && -1 < a && a < 1);
parser.parse(varargin{:});

a1 = parser.Results.InputAlpha;
a2 = parser.Results.OutputAlpha;
if parser.Results.OutputOrder
  m2 = parser.Results.OutputOrder;
else
  m2 = parser.Results.InputOrder;
end

x = eye(parser.Results.InputOrder);
y = srkwii.signal.freqt(x, 'OutputOrder', m2, 'InputAlpha', a1, 'OutputAlpha', a2);
