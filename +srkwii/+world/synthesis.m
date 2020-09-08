function [wav, fs] = synthesis(varargin)
% synthesis - WORLD synthesis
%   wav = synthesis(parameter)
%   wav = synthesis(f0_parameter, sp_parameter)
%   wav = synthesis(fs, temporal_positions, f0, vuv, spectrogram, aperiodicity)
%   wav = synthesis(fs, frame_period, f0, vuv, spectrogram, aperiodicity)
%   wav = synthesis(fs, temporal_positions, f0, vuv, spectrogram, coarse_ap, coarse_axis)
%   wav = synthesis(fs, frame_period, f0, vuv, spectrogram, coarse_ap, coarse_axis)
%   wav = synthesis(___, Name, Value)
%
% Options:
%   DefaultF0: (default: 500)
%   Requiem: use SynthesisRequiem (default: true)
%   RequiemFFTSize: (default: 1024 * 2 ^ ceil(log2(fs / 48000)))
%   RequiemNoiseLength: (default: 2 ^ ceil(log2(fs / 2)))
%   RequiemWindow: (default: (0 : fft_size / 2) * fs / fft_size)
%   RequiemFrequencyRange: (default: (coarse_axis(2) - coarse_axis(1)) * 2, or 6000)
%   ZeroPhase: use zero-phase filtering (available only for Requiem) (default: true)
% Output:
%   wav [T]: synthesized waveform
%   fs: sampling frequency (for convinience)

nargs = 0;
for index = 1:numel(varargin)
  if isstring(varargin{index}) || ischar(varargin{index})
    break; end
  nargs = nargs + 1;
end

parser = inputParser;
if nargs == 1
  addRequired(parser, 'parameter', @isstruct);
  parameter = varargin{1};
  fs = parameter.fs;
  tp = parameter.temporal_positions;
  f0 = parameter.f0;
  if isfield(parameter, 'vuv')
    vu = parameter.vuv;
  else
    vu = [];
  end
  sp = parameter.spectrogram;
  if isfield(parameter, 'coarse_ap') && isfield(parameter, 'coarse_axis')
    bap = parameter.coarse_ap;
    bap_axis = parameter.coarse_axis;
  else
    ap = parameter.aperiodicity;
  end
elseif nargs == 2
  addRequired(parser, 'f0_parameter', @isstruct);
  addRequired(parser, 'sp_parameter', @isstruct);
  f0_parameter = varargin{1};
  sp_parameter = varargin{2};
  fs = sp_parameter.fs;
  tp = f0_parameter.temporal_positions;
  f0 = f0_parameter.f0;
  if isfield(f0_parameter, 'vuv')
    vu = f0_parameter.vuv;
  else
    vu = [];
  end
  sp = sp_parameter.spectrogram;
  if isfield(f0_parameter, 'coarse_ap') && isfield(f0_parameter, 'coarse_axis')
    bap = f0_parameter.coarse_ap;
    bap_axis = f0_parameter.coarse_axis;
  else
    ap = f0_parameter.aperiodicity;
  end
elseif nargs == 6
  addRequired(parser, 'fs', @isscalar);
  addRequired(parser, 'temporal_positions', @isreal);
  addRequired(parser, 'f0', @isreal);
  addRequired(parser, 'vuv', @isreal);
  addRequired(parser, 'spectrogram', @isreal);
  addRequired(parser, 'aperiodicity', @isreal);
  fs = varargin{1};
  tp = varargin{2};
  f0 = varargin{3};
  vu = varargin{4};
  sp = varargin{5};
  ap = varargin{6};
  if numel(tp) == 1
    tp = (0:numel(f0) - 1) * tp / 1000; end
elseif nargs == 7
  addRequired(parser, 'fs', @isscalar);
  addRequired(parser, 'temporal_positions', @isreal);
  addRequired(parser, 'f0', @isreal);
  addRequired(parser, 'vuv', @isreal);
  addRequired(parser, 'spectrogram', @isreal);
  addRequired(parser, 'coarse_ap', @isreal);
  addRequired(parser, 'coarse_axis', @isreal);
  fs = varargin{1};
  tp = varargin{2};
  f0 = varargin{3};
  vu = varargin{4};
  sp = varargin{5};
  bap = varargin{6};
  bap_axis = varargin{7};
  if numel(tp) == 1
    tp = (0:numel(f0) - 1) * tp / 1000; end
else
  error('incorrect usage');
end
addOptional(parser, 'DefaultF0', 500, @isscalar);
addOptional(parser, 'Requiem', true, @isscalar);
addOptional(parser, 'RequiemFFTSize', 1024 * 2 ^ ceil(log2(fs / 48000)), @isscalar);
addOptional(parser, 'RequiemNoiseLength', 2 ^ ceil(log2(fs / 2)), @isscalar);
addOptional(parser, 'RequiemWindow', []);
addOptional(parser, 'RequiemFrequencyRange', []);
addOptional(parser, 'ZeroPhase', true, @isscalar);
parse(parser, varargin{:});

if isempty(vu)
  vu = double(f0 ~= 0);
elseif islogical(vu)
  vu = double(vu);
end

if ~parser.Results.Requiem
  if ~exist('ap', 'var')
    ap = srkwii.world.bap2ap(bap, bap_axis, linspace(0, fs / 2, size(sp, 1))); end
  wav = Synthesis(fs, tp, f0, vu, sp, ap, parser.Results.DefaultF0);
  return
end

% Requiem
fft_size = parser.Results.RequiemFFTSize;
w = parser.Results.RequiemWindow;
if isempty(w)
  w = (0 : fft_size / 2) * fs / fft_size; end
if exist('bap', 'var')
  frequency_interval = bap_axis(2) - bap_axis(1);
  number_of_aperiodicities = numel(bap_axis);
else
  frequency_interval = 3000;
  upper_limit = 15000;
  number_of_aperiodicities = 2 + floor(min(upper_limit, fs / 2 - frequency_interval) / frequency_interval);
  bap_axis = (0:number_of_aperiodicities - 1) * frequency_interval;
  bap = srkwii.world.ap2bap(ap, linspace(0, fs / 2, size(sp, 1)), bap_axis);
end

frequency_range = parser.Results.RequiemFrequencyRange;
if isempty(frequency_range)
  frequency_range = frequency_interval * 2; end

seeds_signals = GetSeedsSignals(...
  fs, fft_size, parser.Results.RequiemNoiseLength, w, frequency_interval, frequency_range, number_of_aperiodicities);
if parser.Results.ZeroPhase
  wav = SynthesisRequiemZero(fs, tp, f0, vu, sp, bap, parser.Results.DefaultF0, seeds_signals);
else
  wav = SynthesisRequiem(fs, tp, f0, vu, sp, bap, parser.Results.DefaultF0, seeds_signals);
end
