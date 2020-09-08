function [fbmatrix, dctmatrix] = mfccmatrix(varargin)
% mfccmatrix - MFCC calculation matrices
%
% Output:
%   fbmatrix [N, L / 2]: filter bank matrix
%                      : [N, L / 2 + 1] if OmitPower is false
%   dctmatrix    [K, N]: DCT matrix
% Options:
%   FFTSize     (=L): (default: 1024)
%   Order       (=K): (default: 12)
%   FilterOrder (=N): (default: 20)
%   SamplingRate    : (default: 16000)
%   OmitPower       : (default: true)

parser = inputParser;
parser.addOptional('FFTSize', 1024, @isscalar);
parser.addOptional('Order', 12, @isscalar);
parser.addOptional('FilterOrder', 20, @isscalar);
parser.addOptional('SamplingRate', 16000, @isscalar);
parser.addOptional('OmitPower', true, @isscalar);
parser.parse(varargin{:});

N = parser.Results.FilterOrder;
K = parser.Results.Order;
L = parser.Results.FFTSize / 2 + 1;
fs = parser.Results.SamplingRate;

mel_axis = hz2mel(linspace(0, fs / 2, L));
interval = mel_axis(end) / (N + 1);
fbmatrix = 1 - abs(mel_axis / interval - (1:N)');
fbmatrix(fbmatrix < 0) = 0;
if parser.Results.OmitPower
  fbmatrix = fbmatrix(:, 2:end); end

dctmatrix = sqrt(2 / N) * cos((pi / N) * ((1:N) - 0.5) .* (1:K)');


function mel = hz2mel(hz)
mel = 1127.01048 * log(hz / 700 + 1);
