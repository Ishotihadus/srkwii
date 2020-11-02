function melsp = sp2melsp(sp, fs, order, reverse)
% sp2melsp - convert spectrogram to mel-spectrogram
%   melsp = sp2melsp(sp, fs, order)
%   melsp = sp2melsp(sp, fs, order, reverse)
%
% Input:
%   sp  [N, T]: log spectrogram
%               mel-scaled log spectrogram if reverse
%   fs        : sampling frequency
%   order (=D): order of mel spectrogram
% Output:
%   melsp [D+1, T]: mel-scaled log spectrogram
%                   linear log spectrogram if reverse
%
% Note: sp should be log scale.

if nargin <= 3
  reverse = false; end

if reverse
  xq = linspace(0, fs, size(sp, 1));
  x = linspace(0, srkwii.signal.mel2hz(xq(end)), order + 1);
else
  x = linspace(0, fs, size(sp, 1));
  xq = linspace(0, srkwii.signal.hz2mel(x(end)), order + 1);
end
melsp = interp1(x, sp, xq, 'linear');
