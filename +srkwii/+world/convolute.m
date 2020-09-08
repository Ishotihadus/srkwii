function y = convolute(x, fs, temporal_positions, spectrogram, varargin)
% convolute - WORLD convolution a.k.a. SP-WORLD
%   y = convolute(x, fs, temporal_positions, spectrogram)
%   y = convolute(x, fs, frame_period, spectrogram)
%   y = convolute(___, Name, Value)
%
% Options:
%   ZeroPhase: (default: false)
%
% References:
%   H. Suda, et al. A revisit to feature handling for high-quality voice conversion based on Gaussian mixture model. In Proc. 2018 APSIPA ASC, pp. 816–822, 2018.

parser = inputParser;
addRequired(parser, 'x', @isreal);
addRequired(parser, 'fs', @isscalar);
addRequired(parser, 'temporal_positions', @isreal);
addRequired(parser, 'spectrogram', @isreal);
addOptional(parser, 'ZeroPhase', false, @isscalar);
parse(parser, x, fs, temporal_positions, spectrogram, varargin{:});

fft_size = (size(spectrogram, 1) - 1) * 2;
spectrogram = log(max(spectrogram, eps));

if parser.Results.ZeroPhase
  spectrogram = spectrogram / 4;
else
  spectrogram = spectrogram / 2;
end

num_samples = numel(x);
y = zeros(num_samples + fft_size, 1);

if numel(temporal_positions) == 1
  temporal_positions = (0:size(spectrogram, 2) - 1) * temporal_positions / 1000; end

dc_remover_base = hanning(fft_size, 'periodic');
dc_remover_base = dc_remover_base / sum(dc_remover_base);

for i = 1 : fft_size
  index = i : fft_size : num_samples;
  time = (index - 1) / fs;
  time = min(max(time, temporal_positions(1)), temporal_positions(end));
  spectrum_slice = interp1(temporal_positions, spectrogram', time, 'linear', 'extrap');
  spectrum_slice = spectrum_slice';
  response = GetPeriodicResponse(spectrum_slice, fft_size);
  dc_remover = dc_remover_base * sum(response, 1);
  response = response - dc_remover;
  buffer_index = (1 : fft_size * numel(index)) + i - 1;
  response = response .* x(index)';
  y(buffer_index) = y(buffer_index) + response(:);
end

y = y(fft_size / 2 + 1 : end - fft_size / 2);

if parser.Results.ZeroPhase
  x = flip(y);
  y = zeros(num_samples + fft_size, 1);
  for i = 1 : fft_size
    index = i : fft_size : num_samples;
    time = (num_samples - index) / fs;
    time = min(max(time, temporal_positions(1)), temporal_positions(end));
    spectrum_slice = interp1(temporal_positions, spectrogram', time, 'linear', 'extrap');
    spectrum_slice = spectrum_slice';
    response = GetPeriodicResponse(spectrum_slice, fft_size);
    dc_remover = dc_remover_base * sum(response, 1);
    response = response - dc_remover;
    buffer_index = (1 : fft_size * numel(index)) + i - 1;
    response = response .* x(index)';
    y(buffer_index) = y(buffer_index) + response(:);
  end
  y = flip(y(fft_size / 2 + 1 : end - fft_size / 2));
end


function response = GetPeriodicResponse(spectrum, fft_size)
spectrum = [spectrum; spectrum(end - 1 : -1 : 2, :)];

cepstrum = 2 * real(fft(spectrum));
cepstrum(1, :) = cepstrum(1, :) / 2;
cepstrum(2 : fft_size / 2, :) = 0;

spectrum = exp(ifft(cepstrum));
spectrum = spectrum(1 : fft_size / 2 + 1, :);
spectrum = [spectrum; conj(spectrum(end - 1 : -1 : 2, :))];
response = fftshift(real(ifft(spectrum)), 1);
