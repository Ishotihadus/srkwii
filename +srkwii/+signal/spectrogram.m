function s = spectrogram(x, window, period)
if nargin < 2
  window = blackman(256); end
window = window(:);
framelen = size(window, 1);

if nargin < 3
  period = 100; end

x = x(:);
x = [zeros(floor(framelen / 2), 1); x; zeros(ceil(framelen / 2), 1)];
len = size(x, 1);

numframes = ceil((len - framelen) / period);
filllen = numframes * period + framelen - len;
x = [x; zeros(filllen, 1)];

s = cell2mat(arrayfun(@(n) x((n * period + 1):(n * period + framelen)), 0:numframes, 'un', 0));
s = fft(s .* window, framelen, 1);

l = floor(framelen / 2) + 1;
s = s(1:l, :);
