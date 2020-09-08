function [min, max] = minmax(input, varargin)
% minmax - find minimum and maximum values (for test)
%   output = minmax(input)
%   output = minmax(___, Name, Value)

parser = inputParser;
parser.addRequired('input');
parser.addOptional('NBest', 1, @(a) isscalar(a) && 1 <= a);
parser.addOptional('Parallel', false, @isscalar);
parser.parse(input, varargin{:});

nbest = parser.Results.NBest;

builder = srkwii.sptk.CommandBuilder('minmax');
builder.AddOption('l', 'int', size(input, 1));
builder.AddOption('b', 'int', nbest);
builder.AddOption('o', 'int', 0);

if parser.Results.Parallel
  output = builder.Exec1by1Parallel(input, size(input, 1) * 4);
else
  output = builder.Exec1by1(input);
end
output = reshape(output, nbest * 2, []);
min = output(1:nbest, :);
max = output(nbest + 1:end, :);
