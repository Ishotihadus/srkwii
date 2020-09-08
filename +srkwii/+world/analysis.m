function parameter = analysis(varargin)
% analysis - WORLD analysis
%   parameter = analysis(x, fs)
%   parameter = analysis(___, Name, Value)
%
% Input:
%   x [T]: waveform
%   fs: sampling frequency
% Options:
%   FramePeriod: フレームシフト (default: basic_frame_period = 1)
%   F0Floor: 仮定される最小基本周波数 (default: 71)
%   F0Ceil: 仮定される最大基本周波数 (default: 800)
%   BasicFramePeriod: Harvest — 分析フレームシフト (default: 1)
%   TargetFs: Harvest — ターゲットサンプリング周波数 (default: 8000)
%   ChannelsInOctave: Harvest — 1 オクターブあたりの帯域フィルタの数 (default: 40)
%   BlurFrames: Harvest — 時間的なぼかしを行うフレーム数（片側） (default: 3)
%   JumpAllowedRange: Harvest — 隣接フレームで許される基本周波数の跳躍 (default: 0.008)
%   VoiceRangeMinimum: Harvest — 有声区間の最短フレーム数 (default: 6)
%   VoiceExtendAllowedRange: Harvest — 隣接する無声フレームにおいて有効候補となる基本周波数候補の範囲 (default: 0.18)
%   VoiceExtendMaxLength: Harvest — 有声区間延長の最大フレーム数 (default: 100)
%   VoiceExtendMinPeriod: Harvest — 有声区間における平均基本周波数の逆数に対する区間の最小フレーム数の比 (default: 2200)
%   UvThreshold: Harvest — 無声区間の最短フレーム数 (default: 9)
%   FilterB, FilterA: Harvest — 基本周波数に対するローパスフィルタの係数 (default: butter(2, 30 / 500))
%   DCCorrection: CheapTrick — パワースペクトルの DC 補正 (default: true)
%   Q0: CheapTrick — TANDEM のスペクトル補償係数 q0 (default: 1 - 2 * q1 = 1.3)
%   Q1: CheapTrick — TANDEM のスペクトル補償係数 q1 (default: -0.15)
%   DefaultF0: CheapTrick, D4C — 無声区間における F0 の計算上の代表値 (default: 500)
%   FFTSize: CheapTrick, D4C — FFT 長 (default: 2 ^ ceil(log2(3 * fs / f0_floor + 1)))
%   D4CThreshold: D4C — LoveTrain（VUV 判定）における閾値（7900 Hz 以下に対する 4000 Hz 以下のパワー比） (default: 0.85)
%   UpperLimit: D4C — バンド非周期性指標の最大周波数 (default: 15000)
%   FrequencyInterval: D4C — バンド非周期性指標で計算される周波数の周期 (default: 3000)
% Output:
%   parameter: parameters
%     parameter.temporal_positions
%     parameter.f0
%     parameter.vuv
%     parameter.f0_candidates
%     parameter.f0_candidate_score
%     parameter.spectrogram
%     parameter.fs
%     parameter.aperiodicity
%     parameter.coarse_ap
%     parameter.coarse_axis
%     parameter.frequency_axis

parser = inputParser;
if ischar(varargin{1}) || isstring(varargin{1})
  addRequired(parser, 'wavfile');
  [x, fs] = audioread(varargin{1});
else
  addRequired(parser, 'x', @isreal);
  addRequired(parser, 'fs', @isscalar);
  [x, fs] = varargin{1:2};
end
addOptional(parser, 'FramePeriod', [], @isscalar);
addOptional(parser, 'F0Floor', 71, @isscalar);
addOptional(parser, 'F0Ceil', 800, @isscalar);
addOptional(parser, 'BasicFramePeriod', 1, @isscalar);
addOptional(parser, 'TargetFs', 8000, @isscalar);
addOptional(parser, 'ChannelsInOctave', 40, @isscalar);
addOptional(parser, 'BlurFrames', 3, @isscalar);
addOptional(parser, 'JumpAllowedRange', 0.008, @isscalar);
addOptional(parser, 'VoiceRangeMinimum', 6, @isscalar);
addOptional(parser, 'VoiceExtendAllowedRange', 0.18, @isscalar);
addOptional(parser, 'VoiceExtendMaxLength', 100, @isscalar);
addOptional(parser, 'VoiceExtendMinPeriod', 2200, @isscalar);
addOptional(parser, 'VoiceExtendMaxSearchSamples', 4, @isscalar);
addOptional(parser, 'UvThreshold', 9, @isscalar);
addOptional(parser, 'FilterB', [0.0078202080334971724, 0.015640416066994345, 0.0078202080334971724], @isreal);
addOptional(parser, 'FilterA', [1.0, -1.7347257688092754, 0.76600660094326412], @isreal);
addOptional(parser, 'DCCorrection', true, @isscalar);
addOptional(parser, 'Q0', [], @isscalar);
addOptional(parser, 'Q1', -0.15, @isscalar);
addOptional(parser, 'DefaultF0', 500, @isscalar);
addOptional(parser, 'FFTSize', []);
addOptional(parser, 'D4CThreshold', 0.85, @isscalar);
addOptional(parser, 'UpperLimit', 15000, @isscalar);
addOptional(parser, 'FrequencyInterval', 3000, @isscalar);
parse(parser, varargin{:});

harvest_option.frame_period = parser.Results.FramePeriod;
harvest_option.f0_floor = parser.Results.F0Floor;
harvest_option.f0_ceil = parser.Results.F0Ceil;
harvest_option.basic_frame_period = parser.Results.BasicFramePeriod;
if isempty(harvest_option.frame_period)
  harvest_option.frame_period = harvest_option.basic_frame_period; end
harvest_option.target_fs = parser.Results.TargetFs;
harvest_option.channels_in_octave = parser.Results.ChannelsInOctave;
harvest_option.blur_frames = parser.Results.BlurFrames;
harvest_option.jump_allowed_range = parser.Results.JumpAllowedRange;
harvest_option.voice_range_minimum = parser.Results.VoiceRangeMinimum;
harvest_option.v_extend_allowed_range = parser.Results.VoiceExtendAllowedRange;
harvest_option.v_extend_max_length = parser.Results.VoiceExtendMaxLength;
harvest_option.v_extend_min_period = parser.Results.VoiceExtendMinPeriod;
harvest_option.v_extend_max_search = parser.Results.VoiceExtendMaxSearchSamples;
harvest_option.uv_threshold = parser.Results.UvThreshold;
harvest_option.filter_b = parser.Results.FilterB;
harvest_option.filter_a = parser.Results.FilterA;
parameter = Harvest(x, fs, harvest_option);

cheaptrick_option.f0_low_limit = parser.Results.F0Floor;
cheaptrick_option.dc_correction = parser.Results.DCCorrection;
cheaptrick_option.q0 = parser.Results.Q0;
cheaptrick_option.q1 = parser.Results.Q1;
if isempty(cheaptrick_option.q0)
  cheaptrick_option.q0 = 1 - 2 * cheaptrick_option.q1; end
cheaptrick_option.default_f0 = parser.Results.DefaultF0;
if ~isempty(parser.Results.FFTSize)
  cheaptrick_option.fft_size = parser.Results.FFTSize; end
sp = CheapTrick(x, fs, parameter, cheaptrick_option);
parameter.spectrogram = sp.spectrogram;
parameter.fs = fs;

d4c_option.threshold = parser.Results.D4CThreshold;
d4c_option.fft_size = parser.Results.FFTSize;
if isempty(d4c_option.fft_size)
  d4c_option.fft_size = (size(parameter.spectrogram, 1) - 1) * 2; end
d4c_option.upper_limit = parser.Results.UpperLimit;
d4c_option.frequency_interval = parser.Results.FrequencyInterval;
parameter = D4C(x, fs, parameter, d4c_option);
