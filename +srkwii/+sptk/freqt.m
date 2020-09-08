function output = freqt(input, varargin)
% freqt - frequency transformation
%   output = freqt(input)
%   output = freqt(___, Name, Value)
%
% Options:
%   OutputOrder: order of warped sequence (default: same as input)
%   InputAlpha : all-pass constant of input sequence
%   OutputAlpha: all-pass constant of output sequence

parser = inputParser;
parser.addRequired('input', @isreal);
parser.addOptional('OutputOrder', size(input, 1) - 1, @(a) isscalar(a) && 0 <= a);
parser.addOptional('InputAlpha', 0, @(a) isscalar(a) && -1 < a && a < 1);
parser.addOptional('OutputAlpha', 0.35, @(a) isscalar(a) && -1 < a && a < 1);
parser.addOptional('Parallel', false, @isscalar);
parser.parse(input, varargin{:});

builder = srkwii.sptk.CommandBuilder('freqt');
builder.AddOption('m', 'int', size(input, 1) - 1);
builder.AddOption('M', 'int', parser.Results.OutputOrder);
builder.AddOption('a', 'float', parser.Results.InputAlpha);
builder.AddOption('A', 'float', parser.Results.OutputAlpha);

if parser.Results.Parallel
  output = builder.Exec1by1Parallel(input, size(input, 1) * 65536);
else
  output = builder.Exec1by1(input);
end

output = reshape(output, parser.Results.OutputOrder + 1, []);
