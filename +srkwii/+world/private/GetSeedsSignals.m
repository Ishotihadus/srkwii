function seeds_signals = GetSeedsSignals(fs, fft_size, noise_length, w, frequency_interval, frequency_range, number_of_aperiodicities)
% Generation of parameters for SynthesisRequiem()
% seeds_signals = GetSeedsSignals(fs)
% seeds_signals = GetSeedsSignals(fs, option)
%
% Input
%   fs     : sampling frequency
%   option : user can cantrol two parameters (but I don't recommend to
%   change them)
%
% Output
%   seeds_signals : generated parameters
%
% Based on 0.2.3 (2018/04/04)

randstream = RandStream('dsfmt19937', 'Seed', 1);

pulse = zeros(fft_size, number_of_aperiodicities);
noise = zeros(noise_length, number_of_aperiodicities);

modified_velvet_noise = GenerateModifiedVelvetNoise(noise_length, fs, randstream);
spec_n = fft(modified_velvet_noise, noise_length);

% Excitation signals in vocal cord vibrations and aperiodic noise were
% generated
for i = 1 : number_of_aperiodicities
  spec = 0.5 + 0.5 * cos(((w - (frequency_interval * (i - 1))) / frequency_range)  * 2 * pi)';
  spec(w > (frequency_interval * i)) = 0;
  spec(w < (frequency_interval * (i - 2))) = 0;
  if i == number_of_aperiodicities
    spec(w > frequency_interval * (i - 1)) = 1; end

  pulse(:, i) = fftshift(real(ifft([spec; spec(end - 1 : -1 : 2)])));
  noise(:, i) = real(ifft(spec_n .* fft(pulse(:, i), noise_length)));
end
pulse(:, 1) = pulse(:, 1) - mean(pulse(:, 1)) .* hanning(fft_size) / mean(hanning(fft_size));
seeds_signals.pulse = pulse;
seeds_signals.noise = noise;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function n = GenerateModifiedVelvetNoise(N, fs, randstream)
base_period = [8, 30, 60];
short_period = 8 * round(base_period * fs / 48000);
n = zeros(N + max(short_period) + 1, 1);

index = 0;
while 1
  v_len = randstream.randi(length(short_period));
  tmp = GenerateShortVelvetNoise(short_period(v_len), randstream);
  n(index + 1 : index + short_period(v_len)) = tmp;
  index = index + short_period(v_len);
  if index >= N
    break; end
end
n = n(1 : N);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function n = GenerateShortVelvetNoise(N, randstream)
n = zeros(N, 1);
td = 4;
range = round(N / td);

safety_rand = 2 * [ones(range / 2, 1); -ones(range / 2, 1)];
for i = 1 : range
  tmp_index = randstream.randi(range);
  tmp = safety_rand(tmp_index);
  safety_rand(tmp_index) = safety_rand(i);
  safety_rand(i) = tmp;
end
n(td * (0 : range - 1)' + randstream.randi(td, range, 1)) = safety_rand;
