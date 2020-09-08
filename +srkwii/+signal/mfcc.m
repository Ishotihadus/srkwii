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

l = size(input, 1);

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
parser.parse(input, varargin{:});

energy_sum = sum(input .^ 2);
energy = log(energy_sum);
energy(energy_sum <= 0) = -1e10;

preemphasis = parser.Results.Preemphasis;
if preemphasis ~= 0
  input = filter([1 -preemphasis], 1, input, -input(1, :) * preemphasis); end

if numel(parser.Results.Window) == l
  input = input .* parser.Results.Window(:);
elseif parser.Results.Window
  w = 0.54 - 0.46 * cos(2 * pi * linspace(0, 1, l));
  w = max([w(1:ceil(l/2)) fliplr(w(1:floor(l/2)))], 0)';
  input = input .* w;
end

L = parser.Results.FFTLength / 2 + 1;
N = parser.Results.FilterOrder;

spec = abs(fft([input; zeros(parser.Results.FFTLength - l, size(input, 2))]));
spec = spec(1:L, :);

mel_list = linspace(0, hz2mel(parser.Results.SamplingRate / 2), N + 2);
filters = cell2mat(arrayfun(@(n) make_filter(mel_list(n:n+2), L, parser.Results.SamplingRate), 1:N, 'un', 0));
filtered = reshape(sum(reshape(spec, L, 1, []) .* filters), N, []);
mfccs = dct(log(max(parser.Results.EnergyFloor, filtered)));

K = parser.Results.OutputOrder;
c0 = mfccs(1, :) * sqrt(2);
output = mfccs(2:K + 1, :);

lifter = parser.Results.Liftering;
if lifter > 0
  output = (1 + lifter / 2 * sin((1:K)' / lifter * pi)) .* output; end


function output = make_filter(points, L, fs)
mel = hz2mel(linspace(0, fs/2, L)');
output = zeros(L, 1);
output(points(1) <= mel & mel <= points(2)) = (mel(points(1) <= mel & mel <= points(2)) - points(1)) / (points(2) - points(1));
output(points(2) <= mel & mel <= points(3)) = (points(3) - mel(points(2) <= mel & mel <= points(3))) / (points(3) - points(2));


function mel = hz2mel(hz)
mel = 1127.01048 * log(hz / 700 + 1);
