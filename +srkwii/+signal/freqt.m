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

mtrx = srkwii.signal.freqtmtrx(...
  'InputOrder', size(input, 1) - 1, 'OutputOrder', parser.Results.OutputOrder,...
  'InputAlpha', parser.Results.InputAlpha, 'OutputAlpha', parser.Results.OutputAlpha);
y = mtrx * input;
