function output = freqtsp(input, varargin)
% freqtsp - frequency transformation of spectrogram
%   output = freqtsp(input)
%   output = freqtsp(___, Name, Value)
%
% Option:
%   InputAlpha : all-pass constant of input sequence
%   OutputAlpha: all-pass constant of output sequence
%   Format     : output format (dB, ln, abs, power) (default: ln)

parser = inputParser;
parser.addRequired('input', @isreal);
parser.addOptional('InputAlpha', 0, @(a) isscalar(a) && -1 < a && a < 1);
parser.addOptional('OutputAlpha', 0.35, @(a) isscalar(a) && -1 < a && a < 1);
parser.addOptional('Format', 'ln');
parser.parse(input, varargin{:});

sp = input;
fmt = validatestring(parser.Results.Format, {'db', 'ln', 'abs', 'power'}, 'freqtsp', 'Format');

switch fmt
case 'db'
  sp = log(10) / 20 * sp;
case 'ln'
  % pass
case 'abs'
  sp = log(max(1e-100, sp));
case 'power'
  sp = log(max(1e-200, sp)) / 2;
end

cep = real(ifft([sp(1:end, :); sp(end-1:-1:2, :)]));
cep = cep(1 : end / 2 + 1, :);
cep = srkwii.signal.freqt(cep, 'InputAlpha', parser.Results.InputAlpha, 'OutputAlpha', parser.Results.OutputAlpha);

output = real(fft([cep(1:end, :); cep(end-1:-1:2, :)]));
output = output(1 : end / 2 + 1, :);

switch fmt
case 'db'
  output = output * 20 / log(10);
case 'ln'
  % pass
case 'abs'
  output = exp(output);
case 'power'
  output = exp(2 * output);
end
