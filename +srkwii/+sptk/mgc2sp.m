function output = mgc2sp(input, varargin)
% mgc2sp - transform mel-generalized cepstrum to spectrum
%   output = mgc2sp(input)
%   output = mgc2sp(___, Name, Value)
%
% Option:
%   Alpha          : alpha α (-1 < α < 1) (default: 0)
%   Gamma          : power parameter γ of mel-generalized cepstrum (-1 <= γ <= 0) (default: 0)
%   NormalizedInput: regard input as normalized cepstrum (default: false)
%   MultipliedInput: regard input as multiplied by γ (default: false)
%   FFTLength      : FFT length (default: 256)
%   OutputFormat   : output format (dB, ln, abs, power, nrad, rad, or deg) (default: ln)
%   Parallel       : (default: false)

parser = inputParser;
parser.addRequired('input', @isreal);
parser.addOptional('Alpha', 0, @(a) isscalar(a) && -1 < a && a < 1);
parser.addOptional('Gamma', 0, @(a) isscalar(a) && -1 <= a && a <= 0);
parser.addOptional('NormalizedInput', false, @isscalar);
parser.addOptional('MultipliedInput', false, @isscalar);
parser.addOptional('FFTLength', 256, @(a) isscalar(a) && a == 2 ^ nextpow2(a));
parser.addOptional('OutputFormat', 'ln')
parser.addOptional('Parallel', false, @isscalar);
parser.parse(input, varargin{:});

numorder = size(input, 1) - 1;

builder = srkwii.sptk.CommandBuilder('mgc2sp');
builder.AddOption('a', 'float', parser.Results.Alpha);
builder.AddOption('g', 'float', parser.Results.Gamma);
builder.AddOption('m', 'int', numorder);
builder.AddOption('n', 'bool', parser.Results.NormalizedInput);
builder.AddOption('u', 'bool', parser.Results.MultipliedInput);
builder.AddOption('l', 'int', parser.Results.FFTLength);

format_idx = selectstr(parser.Results.OutputFormat,...
  {'db', 'ln', 'abs', 'power', 'nrad', 'rad', 'deg'}, 'mgc2sp', 'OutputFormat');
if format_idx >= 4
  builder.AddOption('p', 'bool', true);
  builder.AddOption('o', 'int', format_idx - 4);
else
  builder.AddOption('o', 'int', format_idx);
end

if parser.Results.Parallel
  output = builder.Exec1by1Parallel(input, size(input, 1) * 1024);
else
  output = builder.Exec1by1(input);
end

output = reshape(output, parser.Results.FFTLength / 2 + 1, []);
