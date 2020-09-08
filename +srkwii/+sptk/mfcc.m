function [output, c0, energy] = mfcc(input, varargin)
% mfcc - mel-frequency cepstral analysis
%   output = mfcc(input)
%   output = mfcc(___, Name, Value)
%
% Input:
%   x [L, T]: framed input signal
% Options:
%   Preemphasis : preemphasis coefficient (default: 0.97)/
%   Liftering   : liftering coefficient (default: 22)
%   EnergyFloor : flooring value for calculating log(x) in filterbank analysis (default: 1)
%   FFTLength   : frame length for FFT (default: 2 ^ nextpow2(L + 1))
%   OutputOrder : order of MFCC (default: 12)
%   FilterOrder : order of channel for mel-filter bank (default: 20)
%   SamplingRate: sampling frequency (default: 16000)
%   Window      : Use Hamming windowing (default: true)
%   DFT         : Use DFT (not FFT) for DCT calculation (default: false)
%   Parallel    : (default: false)

parser = inputParser;
parser.addRequired('input', @isreal);
parser.addOptional('Preemphasis', 0.97, @(a) isscalar(a) && isreal(a));
parser.addOptional('Liftering', 22, @(a) isscalar(a) && isreal(a) && 0 <= a);
parser.addOptional('EnergyFloor', 1, @(a) isscalar(a) && isreal(a));
parser.addOptional('FFTLength', 2 ^ nextpow2(size(input, 1) + 1), @(a) isscalar(a) && a >= 2 && nextpow2(a) == a);
parser.addOptional('OutputOrder', 12, @(a) isscalar(a) && isreal(a) && 0 < a);
parser.addOptional('FilterOrder', 20, @(a) isscalar(a) && isreal(a) && 0 < a);
parser.addOptional('SamplingRate', 16000, @(a) isscalar(a) && isreal(a) && 0 < a);
parser.addOptional('Window', true, @(a) isscalar(a));
parser.addOptional('DFT', false, @(a) isscalar(a));
parser.addOptional('Parallel', false, @isscalar);
parser.parse(input, varargin{:});

if isempty(input)
  output = [];
  return
end

framelen = size(input, 1);
outputorder = parser.Results.OutputOrder;

builder = srkwii.sptk.CommandBuilder('mfcc');
builder.AddOption('a', 'float', parser.Results.Preemphasis);
builder.AddOption('c', 'int', parser.Results.Liftering);
builder.AddOption('e', 'float', parser.Results.EnergyFloor);
builder.AddOption('s', 'float', parser.Results.SamplingRate / 1000);
builder.AddOption('l', 'int', framelen);
builder.AddOption('L', 'int', parser.Results.FFTLength);
builder.AddOption('m', 'int', outputorder);
builder.AddOption('n', 'int', parser.Results.FilterOrder);
if parser.Results.Window
  builder.AddOption('w', 'int', 0);
else
  builder.AddOption('w', 'int', 1);
end
builder.AddOption('d', 'bool', parser.Results.DFT);
builder.AddOption('E', 'bool', nargout >= 3);
builder.AddOption('0', 'bool', nargout >= 2);

if parser.Results.Parallel
  output = builder.Exec1by1Parallel(input, size(input, 1) * 1024);
else
  output = builder.Exec1by1(input);
end

if nargout >= 3
  output = reshape(output, outputorder + 2, []);
  c0 = output(end-1, :);
  energy = output(end, :);
  output = output(1:end-2, :);
elseif nargout >= 2
  output = reshape(output, outputorder + 1, []);
  c0 = output(end, :);
  output = output(1:end-1, :);
else
  output = reshape(output, outputorder, []);
end
