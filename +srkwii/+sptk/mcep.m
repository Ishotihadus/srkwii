function output = mcep(input, varargin)
% mcep - mel cepstral analysis
%   output = mcep(input)
%   output = mcep(___, Name, Value)

parser = inputParser;
parser.addRequired('input', @isreal);
parser.addOptional('Alpha', 0.35, @(a) isscalar(a) && isreal(a));
parser.addOptional('Order', 25, @(a) isscalar(a) && 0 <= a);
parser.addOptional('Format', 'signal');
parser.addOptional('MinimumIteration', 2, @(a) isscalar(a) && 0 <= a);
parser.addOptional('MaximumIteration', 30, @(a) isscalar(a) && 0 <= a);
parser.addOptional('EndCondition', 1e-3, @(a) isscalar(a) && 0 <= a);
parser.addOptional('PeriodogramFlooring', 0, @(a) isscalar(a) && 0 <= a);
parser.addOptional('MinimumDeterminant', 1e-6, @(a) isscalar(a) && 0 <= a);
parser.addOptional('Parallel', false, @isscalar);
parser.parse(input, varargin{:});

if isempty(input)
  output = [];
  return
end

framelen = size(input, 1);
if parser.Results.Format ~= 0
  framelen = (size(input, 1) - 1) * 2; end

builder = srkwii.sptk.CommandBuilder('mcep');
builder.AddOption('a', 'float', parser.Results.Alpha);
builder.AddOption('m', 'int', parser.Results.Order);
builder.AddOption('l', 'int', framelen);
builder.AddOption('q', 'int', selectstr(parser.Results.Format, {'signal', 'db', 'ln', 'abs', 'power'}, 'mcep', 'Format'));
builder.AddOption('i', 'int', parser.Results.MinimumIteration);
builder.AddOption('j', 'int', parser.Results.MaximumIteration);
builder.AddOption('d', 'float', parser.Results.EndCondition);
builder.AddOption('e', 'float', parser.Results.PeriodogramFlooring);
builder.AddOption('f', 'float', parser.Results.MinimumDeterminant);

if parser.Results.Parallel
  output = builder.Exec1by1Parallel(input, size(input, 1) * 1024);
else
  output = builder.Exec1by1(input);
end
output = reshape(output, parser.Results.Order + 1, []);
