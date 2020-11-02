function [h, u, x, div, history] = nmf(y, varargin)
% nmf - Non-negative matrix factorization
%
% Input:
%   y [K, T]: matrix to be factorized
% Output:
%   h [K, N]: estimated base matrix
%   u [N, T]: estimated activation matrix
%   x [K, T]: reconstructed matrix (= h * u)
%   div     : divergence
% Options:
%   BaseNum           (=N) : number of bases. ignored if Base is specified. (default: 100)
%   Base            [N, K] : initial base
%   Activation      [K, T] : initial activation
%   ConstantBase           : not updated elements for bases (default: false)
%                              apply for all elements if true or false
%                              apply for each bases if [1, N]
%                              apply for each rows if [K, 1]
%   ConstantActivation     : not updated elements for activation in the same way as ConstantBase (default: false)
%   ConstantBaseAvail      : whether activate ConstantBase (default: @(~) true)
%   ConstantActivationAvail: whether activate ConstantActivation (default: @(~) true)
%   Iteration              : number of iterations
%   NormalizeBase          : p of generalized vector norm. disabled if falsy
%   Divergence             : divergence function (default: kl)
%                              use predefined divergence if string
%                                (Euclidean, kl = KullbackLeibler, SparseKL, is = ItakuraSaito, LogEuclidean)
%                              use specified function if function (arguments: y, x / returns: element-wise divergence)
%                                requires UpdateBase and UpdateActivation if custom divergence
%   SparseKLBaseAlpha      : parameter for SparseKL divergence
%   SparseKLActivationAlpha: parameter for SparseKL divergence
%   SparseKLWeight         : parameter for SparseKL divergence
%   AdditionalDivergence   : additional divergence for logging. not used for training
%   Sparsity               : criterion of sparsity penalty for activation
%                              use predefined sparsity criterion if string
%                                (l1 = L1Norm, L12QuasiNorm)
%                              use specified function and criterion if struct
%                                struct('criterion', @(u) ~, 'update', struct('order', { 0; -1 }, 'func', { @(u) ~; @(u) ~ }))
%   SparsityWeight         : weight of sparsity penalty (a.k.a. lambda)
%   Density                : criterion of density penalty for base
%                              use predefined sparsity criterion if string
%                                (l2 = L2Norm)
%                              use specified function and criterion if struct
%                                struct('criterion', @(u) ~, 'update', struct('order', { 0; -1 }, 'func', { @(h) ~; @(h) ~ }))
%   DensityWeight          : weight of density penalty
%   ActivityIndropout      : probability of addition of small random value for zero elements in activation in each iteration
%   CepstralRegularization : struct for cepstral regularization
%                              the structure must be
%                                struct('weight', @weight_function, 'mu', [N, K], 'fbmatrix', [A, K], 'dctmatrix', [B, A])
%   UpdateBase             : update function for base. sparsity is ignored if specified.
%   UpdateActivation       : update function for activation. sparsity is ignored if specified.
%   UseGPU                 : whether to use gpu (default: isa(y, 'gpuArray'))
%   PrintLog               : whether to print log. print per n iteration if number specified

parser = inputParser;
addParameter(parser, 'BaseNum', 100);
addParameter(parser, 'Base', []);
addParameter(parser, 'Activation', []);
addParameter(parser, 'ConstantBase', false);
addParameter(parser, 'ConstantActivation', false);
addParameter(parser, 'ConstantBaseAvail', @(~) true);
addParameter(parser, 'ConstantActivationAvail', @(~) true);
addParameter(parser, 'Iteration', 100, @isscalar);
addParameter(parser, 'NormalizeBase', 1, @isscalar);
addParameter(parser, 'CepstralNormalize', false, @isscalar);
addParameter(parser, 'Divergence', 'kl');
addParameter(parser, 'SparseKLBaseAlpha', 0);
addParameter(parser, 'SparseKLActivationAlpha', 0);
addParameter(parser, 'SparseKLWeight', 1);
addParameter(parser, 'AdditionalDivergence', []);
addParameter(parser, 'Sparsity', 'none');
addParameter(parser, 'SparsityWeight', 1);
addParameter(parser, 'Density', 'none');
addParameter(parser, 'DensityWeight', 1);
addParameter(parser, 'ActivityIndropout', 0);
addParameter(parser, 'CepstralRegularization', []);
addParameter(parser, 'UpdateBase', []);
addParameter(parser, 'UpdateActivation', []);
addParameter(parser, 'UseGPU', isa(y, 'gpuArray'), @isscalar);
addParameter(parser, 'PrintLog', true, @isscalar);
parse(parser, varargin{:});

% データ格納用構造体
d = srkwii.util.HandleStruct('y', y, 'K', size(y, 1), 'T', size(y, 2));

% 乱数生成
randstream = RandStream('dsfmt19937');

% ログ出力の設定
print_log = parser.Results.PrintLog * 1;

% 基底の正規化の p
base_normalize = parser.Results.NormalizeBase;
if ~base_normalize
  base_normalize = false;
else
  base_normalize = base_normalize * 1;
end

% 初期値の生成と大きさの決定
nmf_init(d, parser, randstream);

% 固定成分の決定
h_upd_elem = true;
if numel(parser.Results.ConstantBase) == 1
  if parser.Results.ConstantBase
    h_upd_elem = []; end
else
  h_upd_elem = true(size(parser.Results.ConstantBase));
  h_upd_elem(parser.Results.ConstantBase) = false;
  if ~any(h_upd_elem)
    h_upd_elem = []; end
end
u_upd_elem = true;
if numel(parser.Results.ConstantActivation) == 1
  if parser.Results.ConstantActivation
    u_upd_elem = []; end
else
  u_upd_elem = true(size(parser.Results.ConstantActivation));
  u_upd_elem(parser.Results.ConstantActivation) = false;
  if ~any(u_upd_elem)
    u_upd_elem = []; end
end

% GPU
use_gpu = parser.Results.UseGPU;
if use_gpu
  if ~isa(d.y, 'gpuArray')
    d.y = gpuArray(d.y); end
  if ~isa(d.h, 'gpuArray')
    d.h = gpuArray(d.h); end
  if ~isa(d.u, 'gpuArray')
    d.u = gpuArray(d.u); end
else
  if isa(d.y, 'gpuArray')
    d.y = gather(d.y); end
  if isa(d.h, 'gpuArray')
    d.h = gather(d.h); end
  if isa(d.u, 'gpuArray')
    d.u = gather(d.u); end
end

% ダイバージェンス関数と更新式の決定
fun_div = parser.Results.Divergence;
fun_sparse = [];
fun_dense = [];
fun_cep = [];
if isa(fun_div, 'function_handle')
  % 直接指定されたときはそれに従う
  update_h = parser.Results.UpdateBase;
  update_u = parser.Results.UpdateActivation;
else
  % そうでないときには生成する
  heq = eqs_init;
  ueq = eqs_init;
  if isnumeric(fun_div)
    b = fun_div;
    a1 = b - 1;
    a2 = b - 2;

    if b == 0
      % Itakura-Saito
      fun_div = @(y, x) y ./ x - log(y) + log(x) - 1;
    elseif b == 1
      % Kullback-Leibler
      fun_div = @(y, x) y .* (log(y) - log(x)) - y + x;
    elseif b == 2
      % Euclidean
      fun_div = @(y, x) (y - x).^2 / 2;
    else
      fun_div = @(y, x) (y.^b - y .* x.^a1) / a1 - (y.^b - x.^b) / b;
    end

    if b < 1
      heq = eqs_plus(heq, 0, @(y, x, h, u) x.^a1 * u');
      ueq = eqs_plus(ueq, 0, @(y, x, h, u) h' * x.^a1);
    else
      heq = eqs_plus(heq, a1, @(y, x, h, u) x.^a1 * u' ./ h.^a1);
      ueq = eqs_plus(ueq, a1, @(y, x, h, u) h' * x.^a1 ./ u.^a1);
    end

    if b < 2
      heq = eqs_minus(heq, a2, @(y, x, h, u) y .* x.^a2 * u' ./ h.^a2);
      ueq = eqs_minus(ueq, a2, @(y, x, h, u) h' * (y .* x.^a2) ./ u.^a2);
    else
      heq = eqs_minus(heq, 0, @(y, x, h, u) y .* x.^a2 * u');
      ueq = eqs_minus(ueq, 0, @(y, x, h, u) h' * (y .* x.^a2));
    end
  else
    switch validatestring(fun_div, {'Euclidean', 'kl', 'KullbackLeibler', 'is', 'ItakuraSaito', 'LogEuclidean', 'SparseKL'})
    case 'Euclidean'
      fun_div = @(y, x) (y - x).^2;
      heq = eqs_minus(heq, 0, @(y, x, h, u) 2 * y * u');
      heq = eqs_plus(heq, 1, @(y, x, h, u) 2 * x * u' ./ h);
      ueq = eqs_minus(ueq, 0, @(y, x, h, u) 2 * h' * y);
      ueq = eqs_plus(ueq, 1, @(y, x, h, u) 2 * h' * x ./ u);
    case {'kl', 'KullbackLeibler'}
      fun_div = @(y, x) y .* (log(y) - log(x)) - y + x;
      heq = eqs_minus(heq, -1, @(y, x, h, u) y ./ x * u' .* h);
      heq = eqs_plus(heq, 0, @(y, x, h, u) sum(u, 2)');
      ueq = eqs_minus(ueq, -1, @(y, x, h, u) h' * (y ./ x) .* u);
      ueq = eqs_plus(ueq, 0, @(y, x, h, u) sum(h, 1)');
    case 'SparseKL'
      d.skl_ah = parser.Results.SparseKLBaseAlpha;
      d.skl_au = parser.Results.SparseKLActivationAlpha;
      d.skl_w = parser.Results.SparseKLWeight;
      if ~isa(d.skl_ah, 'function_handle')
        skl_ah = d.skl_ah;
        d.skl_ah = @(~) skl_ah;
      end
      if ~isa(d.skl_au, 'function_handle')
        skl_au = d.skl_au;
        d.skl_au = @(~) skl_au;
      end
      if ~isa(d.skl_w, 'function_handle')
        skl_w = d.skl_w;
        d.skl_w = @(~) skl_w;
      end
      fun_div = @(y, x) y .* (log(y) - log(x)) - y + x;
      heq = eqs_minus(heq, -1, @(y, x, h, u, e) ((y ./ x * u').^e.skl_w(e) .* h).^(1 + e.skl_ah(e)));
      heq = eqs_plus(heq, 0, @(y, x, h, u, e) (sum(u, 2)').^(1 + e.skl_ah(e)));
      ueq = eqs_minus(ueq, -1, @(y, x, h, u, e) ((h' * (y ./ x)).^e.skl_w(e) .* u).^(1 + e.skl_au(e)));
      ueq = eqs_plus(ueq, 0, @(y, x, h, u, e) (sum(h, 1)').^(1 + e.skl_au(e)));
    case {'is', 'ItakuraSaito'}
      fun_div = @(y, x) y ./ x - log(y) + log(x) - 1;
      heq = eqs_plus(heq, 0, @(y, x, h, u) 1 ./ x * u');
      heq = eqs_minus(heq, -2, @(y, x, h, u) y ./ x.^2 * u' .* h.^2);
      ueq = eqs_plus(ueq, 0, @(y, x, h, u) h' * (1 ./ x));
      ueq = eqs_minus(ueq, -2, @(y, x, h, u) h' * (y ./ x.^2) .* u.^2);
    case 'LogEuclidean'
      fun_div = @(y, x) (log(y) - log(x)).^2;
      heq = eqs_plus(heq, 0, @(y, x, h, u) (2 * log(x) + 1 ./ x - 2 * (y < 1) .* log(y)) ./ x * u');
      heq = eqs_minus(heq, -1, @(y, x, h, u) 2 * (y > 1) .* log(y) ./ x * u' .* h);
      heq = eqs_minus(heq, -2, @(y, x, h, u) 1 ./ x.^2 * u' .* h.^2);
      ueq = eqs_plus(ueq, 0, @(y, x, h, u) h' * ((2 * log(x) + 1 ./ x - 2 * (y < 1) .* log(y)) ./ x));
      ueq = eqs_minus(ueq, -1, @(y, x, h, u) 2 * h' * ((y > 1) .* log(y) ./ x) .* u);
      ueq = eqs_minus(ueq, -2, @(y, x, h, u) h' * (1 ./ x.^2) .* u.^2);
    end
  end

  % sparsity penalties
  % 重み関数の決定
  sparse_weight = parser.Results.SparsityWeight;
  if ~isa(sparse_weight, 'function_handle')
    sparse_weight = @(~) sparse_weight; end
  if isstruct(parser.Results.Sparsity)
    % 微分の形式で与えられていればそれに従う
    fun_sparse = parser.Results.Sparsity.criterion;
    s = parser.Results.Sparsity.update;
    for n = 1:numel(s)
      ueq = eqs_plus(ueq, s(n).order, @(y, x, h, u, e) sparse_weight(e) * s(n).func(y, x, h, u, e)); end
  else
    % そうでなければ規定のスパースネスを使う
    switch validatestring(parser.Results.Sparsity, {'none', 'l1', 'L1Norm', 'L12QuasiNorm'})
    case {'l1', 'L1Norm'}
      % \sum_{k,t} u -> 1
      fun_sparse = @(u) sum(u, 1);
      ueq = eqs_plus(ueq, 0, @(~, ~, ~, ~, e) sparse_weight(e));
    case 'L12QuasiNorm' % L1/2 quasi-norm
      % \sum_t (\sum_k \sqrt u)^2 <= \sum_{k,t} u / l (l = sqrt(u) / \sum_k sqrt(u))
      fun_sparse = @(u) sum(sqrt(u), 1) .^ 2;
      ueq = eqs_plus(ueq, 0, @(~, ~, ~, u, e) sum(sqrt(u), 1) ./ (sqrt(u) + eps) .* sparse_weight(e));
    end
  end

  % density penalties
  density_weight = parser.Results.DensityWeight;
  if ~isa(density_weight, 'function_handle')
    density_weight = @(~) density_weight; end
  if isstruct(parser.Results.Density)
    % 微分の形式で与えられていればそれに従う
    fun_dense = parser.Results.Density.criterion;
    s = parser.Results.Density.update;
    for n = 1:numel(s)
      heq = eqs_plus(heq, s(n).order, @(y, x, h, u, e) density_weight(e) * s(n).func(y, x, h, u, e)); end
  else
    % そうでなければ規定のスパースネスを使う
    switch validatestring(parser.Results.Density, {'none', 'L2Norm'})
    case 'L2Norm'
      % \sum_{k,t} h^2 -> 2 h
      fun_dense = @(h) sum(h .^ 2, 1);
      heq = eqs_plus(heq, 1, @(~, ~, h, ~, e) 2 * density_weight(e) * h);
    end
  end

  cep_opts = parser.Results.CepstralRegularization;
  if ~isempty(cep_opts)
    cep_weight = @(~) 0;
    if isfield(cep_opts, 'weight')
      cep_weight = cep_opts.weight; end
    if ~isa(cep_weight, 'function_handle')
      cep_weight = @(~) cep_weight; end
    d.cep_f = cep_opts.fbmatrix;
    d.cep_c = cep_opts.dctmatrix;
    d.cep_mu = cep_opts.mu;
    d.cep_M = size(d.cep_c, 1);
    d.cep_N = size(d.cep_c, 2); % = size(d.cep_f, 1);
    if use_gpu
      if ~isa(d.cep_mu, 'gpuArray')
        d.cep_mu = gpuArray(d.cep_mu); end
      if ~isa(d.cep_f, 'gpuArray')
        d.cep_f = gpuArray(d.cep_f); end
      if ~isa(d.cep_c, 'gpuArray')
        d.cep_c = gpuArray(d.cep_c); end
    else
      if isa(d.cep_mu, 'gpuArray')
        d.cep_mu = gather(d.cep_mu); end
      if isa(d.cep_f, 'gpuArray')
        d.cep_f = gather(d.cep_f); end
      if isa(d.cep_c, 'gpuArray')
        d.cep_c = gather(d.cep_c); end
    end
    cepstrum_first(d);

    fun_cep = @(e) sum((e.cep_c * log(e.cep_f * e.h) - e.cep_mu).^2, 1);
    % -ƒ' * (A / G.^2) .* H.^2
    heq = eqs_minus(heq, -2, @(~, ~, h, ~, e) cep_weight(e) * e.cep_f' * (e.cep_a ./ e.cep_g.^2) .* h.^2);
    % f' * (B < 0 -> B) ./ G .* H
    heq = eqs_plus(heq, -1, @(~, ~, h, ~, e) cep_weight(e) * e.cep_f' * ((e.cep_b < 0) .* e.cep_b ./ e.cep_g) .* h);
    %  f' * (A .* (2logG + 1./G) + (B > 0 -> B)) ./ G
    heq = eqs_plus(heq, 0, @(~, ~, h, ~, e) cep_weight(e) * e.cep_f' * ((e.cep_a .* (2 * log(e.cep_g) + 1 ./ e.cep_g) + (e.cep_b > 0) .* e.cep_b) ./ e.cep_g));
  end

  [update_h, update_h_str] = eqs_compile(heq, 5);
  [update_u, update_u_str] = eqs_compile(ueq, 5);
  if print_log > 0
    fprintf('H update equation:\n');
    eqs_print(heq);
    fprintf('  -> %s\n', update_h_str);
    fprintf('U update equation:\n');
    eqs_print(ueq);
    fprintf('  -> %s\n', update_u_str);
  end
end

% 追加ダイバージェンス
fun_div2 = parser.Results.AdditionalDivergence;
if ~isempty(fun_div2) && ~isa(fun_div2, 'function_handle')
  if isnumeric(fun_div2)
    b = fun_div2;
    if b == 0
      % Itakura-Saito
      fun_div2 = @(y, x) y ./ x - log(y) + log(x) - 1;
    elseif b == 1
      % Kullback-Leibler
      fun_div2 = @(y, x) y .* (log(y) - log(x)) - y + x;
    elseif b == 2
      % Euclidean
      fun_div2 = @(y, x) (y - x).^2 / 2;
    else
      a1 = b - 1;
      fun_div2 = @(y, x) (y.^b - y .* x.^a1) / a1 - (y.^b - x.^b) / b;
    end
  else
    switch validatestring(fun_div2,...
      {'Euclidean', 'kl', 'KullbackLeibler', 'rkl', 'ReverseKullbackLeibler', 'is', 'ItakuraSaito', 'LogEuclidean'})
    case 'Euclidean'
      fun_div2 = @(y, x) (y - x).^2;
    case {'kl', 'KullbackLeibler'}
      fun_div2 = @(y, x) y .* log(y ./ x) - y + x;
    case {'rkl', 'ReverseKullbackLeibler'}
      fun_div2 = @(x, y) y .* log(y ./ x) - y + x;
    case {'is', 'ItakuraSaito'}
      fun_div2 = @(y, x) y ./ x - log(y ./ x) - 1;
    case 'LogEuclidean'
      fun_div2 = @(y, x) (log(y) - log(x)).^2;
    end
  end
end

% history save
save_history = false;
if nargout >= 5
  save_history = true;
  fields = {'iter' 'dist'};
  if isa(fun_div2, 'function_handle')
    fields(end + 1) = 'dist2'; end
  if isa(fun_sparse, 'function_handle')
    fields(end + 1) = 'sparsity'; end
  if isa(fun_dense, 'function_handle')
    fields(end + 1) = 'density'; end
  fields(2, :) = arrayfun(@(~) {}, 1:numel(fields), 'un', 0);
  history = struct(fields{:});
end

% print header
if print_log > 0
  fprintf(' iter   divergence');
  if isa(fun_div2, 'function_handle')
    fprintf('  divergence2'); end
  if isa(fun_sparse, 'function_handle')
    fprintf('     sparsity'); end
  if isa(fun_dense, 'function_handle')
    fprintf('     density'); end
  if ~isempty(h_upd_elem)
    fprintf('  |   divergence');
    if isa(fun_div2, 'function_handle')
      fprintf('  divergence2'); end
    if isa(fun_sparse, 'function_handle')
      fprintf('     sparsity'); end
    if isa(fun_dense, 'function_handle')
      fprintf('     density'); end
    if isa(fun_cep, 'function_handle')
      fprintf(' cepstral div'); end
  end
  fprintf('\n');
end

if parser.Results.CepstralNormalize
  d.y_cep = d.y;
  logdy = log(d.y);
  c0 = mean(gather(logdy), 1);
  d.y = exp(logdy - c0);
  d.u = d.u ./ exp(c0);
end


d.x = d.h * d.u;
d.iter = 0;

calc_div(d, fun_div, fun_div2, fun_sparse, fun_dense, fun_cep);
if print_log > 0
  fprintf(' init: %s\n', diststr(d)); end
if save_history
  history(end + 1) = copy_hist(d); end

for iter = 1:parser.Results.Iteration
  d.iter = iter;
  iter_print_log = mod(iter, print_log) == 0;

  calc_div(d, fun_div, fun_div2, fun_sparse, fun_dense, []);
  if iter_print_log
    fprintf('%5d: %s', iter, diststr(d, true)); end
  if save_history
    history(end + 1) = copy_hist(d); end

  % 基底の更新
  force_update_base = ~parser.Results.ConstantBaseAvail(d);
  if ~isempty(h_upd_elem) || force_update_base
    if ~isempty(cep_opts)
      cepstrum_before(d); end

    if numel(h_upd_elem) == 1 || force_update_base
      d.h = update_h(d.y, d.x, d.h, d.u, d);
    elseif size(h_upd_elem, 2) == 1
      d.h(h_upd_elem, :) = update_h(d.y(h_upd_elem, :), d.x(h_upd_elem, :), d.h(h_upd_elem, :), d.u, d);
    elseif size(h_upd_elem, 1) == 1
      d.h(:, h_upd_elem) = update_h(d.y, d.x, d.h(:, h_upd_elem), d.u(h_upd_elem, :), d);
    else
      h_tmp = update_h(d.y, d.x, d.h, d.u, d);
      d.h(h_upd_elem) = h_tmp(h_upd_elem);
    end

    % 基底のノーマライズ
    if base_normalize || force_update_base
      if numel(h_upd_elem) == 1 || force_update_base
        hnorm = vecnorm(d.h, base_normalize);
        d.h = d.h ./ hnorm;
        d.u = d.u .* hnorm';
      elseif size(h_upd_elem, 1) == 1
        hnorm = vecnorm(d.h(:, h_upd_elem), base_normalize);
        d.h(:, h_upd_elem) = d.h(:, h_upd_elem) ./ hnorm;
        d.u(h_upd_elem, :) = d.u(h_upd_elem, :) .* hnorm';
      end
    end

    d.x = d.h * d.u;

    if ~isempty(cep_opts)
      cepstrum_after(d); end

    calc_div(d, fun_div, fun_div2, fun_sparse, fun_dense, fun_cep);
    if iter_print_log
      fprintf('  |  %s', diststr(d)); end
    if save_history
      history(end + 1) = copy_hist(d); end
  end

  % 生起状態の更新
  force_update_act = ~parser.Results.ConstantActivationAvail(d);
  if ~isempty(u_upd_elem) || force_update_act
    if numel(u_upd_elem) == 1 || force_update_act
      d.u = update_u(d.y, d.x, d.h ,d.u, d);
      if parser.Results.ActivityIndropout > 0
        zero_idx = find(d.u == 0);
        zero_idx = randsample(zero_idx, round(numel(zero_idx) * parser.Results.ActivityIndropout));
        if numel(zero_idx) > 0
          [~, zero_idx_c] = ind2sub(size(d.u), zero_idx);
          val = mean(d.u, 1)';
          d.u(zero_idx) = val(zero_idx_c) .* randstream.rand(size(zero_idx_c));
        end
      end
    elseif size(u_upd_elem, 2) == 1
      d.u(u_upd_elem, :) = update_u(d.y, d.x, d.h(:, u_upd_elem) ,d.u(u_upd_elem, :), d);
    elseif size(u_upd_elem, 1) == 1
      d.u(:, u_upd_elem) = update_u(d.y(:, u_upd_elem), d.x(:, u_upd_elem), d.h ,d.u(:, u_upd_elem), d);
    else
      u_tmp = update_u(d.y, d.x, d.h ,d.u, d);
      d.u(u_upd_elem) = u_tmp(u_upd_elem);
    end

    d.x = d.h * d.u;
  end

  if iter_print_log
    fprintf('\n'); end
end

d.iter = iter + 1;
calc_div(d, fun_div, fun_div2, fun_sparse, fun_dense, fun_cep);
if iter_print_log
  fprintf('final: %s\n', diststr(d, true)); end
if save_history
  history(end + 1) = copy_hist(d); end

if parser.Results.CepstralNormalize
  logdx = log(d.x);
  c0 = c0 - mean(gather(logdx), 1);
  d.x = exp(logdx + c0);
  d.u = d.u .* exp(c0);
  d.y = d.y_cep;

  if iter_print_log
    fprintf('cepfin:%s\n', diststr(d, true)); end
end

h = d.h;
u = d.u;
x = d.x;
div = d.dist;



function nmf_init(d, parser, randstream)
d.h = parser.Results.Base;
d.u = parser.Results.Activation;

if isempty(d.h)
  if isempty(d.u)
    d.N = parser.Results.BaseNum;
  else
    d.N = size(d.u, 1);
  end
else
  d.N = size(d.h, 2);
end

if isempty(d.h)
  d.h = randstream.rand(d.K, d.N); end
if isempty(d.u)
  d.u = randstream.rand(d.N, d.T); end


function cepstrum_first(d)
% create b
if isa(d.y, 'gpuArray')
  d.cep_beta = abs(randn(d.cep_M, d.cep_N, d.N, 'gpuArray'));
else
  d.cep_beta = abs(randn(d.cep_M, d.cep_N, d.N));
end
d.cep_beta = d.cep_beta ./ sum(d.cep_beta, 2);
% G [n, k] = f * H
d.cep_g = d.cep_f * d.h;

function cepstrum_before(d)
% A [n, k] = sum(c.^2 ./ beta, m->)
d.cep_a = reshape(sum(d.cep_c.^2 ./ d.cep_beta, 1), d.cep_N, d.N);
% B [n, k] = sum(c .* gamma ./ beta, m->)
%          = sum(c .^ 2 .* logG ./ beta + [m, n] c .* ([m, k] mu - c * logG))
%          = [n, k] logG .* [n, k] A + [n, m] c' * ([m, k] mu - c * logG)
d.cep_b = log(d.cep_g) .* d.cep_a + d.cep_c' * (d.cep_mu - d.cep_c * log(d.cep_g));

function cepstrum_after(d)
% G [n, k] = f * H
cep_g = d.cep_f * d.h;
% beta = abs(c .* (logG - logG-) - beta .* (mu - c * logG-))
d.cep_beta = abs(...
  d.cep_c .* reshape(log(cep_g) - log(d.cep_g), 1, d.cep_N, d.N)...
  - d.cep_beta .* reshape(d.cep_mu - d.cep_c * log(d.cep_g), d.cep_M, 1, d.N));
d.cep_beta = d.cep_beta ./ sum(d.cep_beta, 2);
d.cep_g = cep_g;
% [m, k] = mu - c * logG
d.cep_mu_minus_clogg = d.cep_mu - d.cep_c * log(cep_g);


% 現在のダイバージェンスを計算する
function calc_div(d, fun_div, fun_div2, fun_sparse, fun_dense, fun_cep)
d.dist = mean(fun_div(d.y, d.x), 'all');
if isa(d.dist, 'gpuArray')
  d.dist = gather(d.dist); end
if isa(fun_div2, 'function_handle')
  d.dist2 = mean(fun_div2(d.y, d.x), 'all');
  if isa(d.dist2, 'gpuArray')
    d.dist2 = gather(d.dist2); end
end
if isa(fun_sparse, 'function_handle')
  d.sparsity = mean(fun_sparse(d.u)) ./ d.T;
  if isa(d.sparsity, 'gpuArray')
    d.sparsity = gather(d.sparsity); end
end
if isa(fun_dense, 'function_handle')
  d.density = mean(fun_dense(d.h)) ./ d.N;
  if isa(d.density, 'gpuArray')
    d.density = gather(d.density); end
end
if isa(fun_cep, 'function_handle')
  d.cep_dist = mean(fun_cep(d));
  if isa(d.cep_dist, 'gpuArray')
    d.cep_dist = gather(d.cep_dist); end
end

% ログ出力用の文字列を生成する
function str = diststr(d, skip_cep)
str = sprintf('%#1.5e', d.dist);
if isfield(d, 'dist2')
  str = sprintf('%s  %#1.5e', str, d.dist2); end
if isfield(d, 'sparsity')
  str = sprintf('%s  %#1.5e', str, d.sparsity); end
if (nargin <= 1 || ~skip_cep) && isfield(d, 'cep_dist')
  str = sprintf('%s  %#1.5e', str, d.cep_dist); end

% ヒストリーを返す用の構造体を作る
function s = copy_hist(d)
s = struct('iter', d.iter, 'dist', d.dist);
if isfield(d, 'dist2')
  s.dist2 = d.dist2; end
if isfield(d, 'sparsity')
  s.sparsity = d.sparsity; end
if isfield(d, 'density')
  s.density = d.density; end
if isfield(d, 'cep_dist')
  s.cep_dist = d.cep_dist; end


% 更新式方程式の初期化
function equations = eqs_init
equations = struct('order', {}, 'func', {}, 'minus', {});

% 更新式方程式の特定次数に項を追加
function equations = eqs_plus(equations, order, func)
orders = arrayfun(@(s) s.order, equations);
idx = find(orders == order, 1);
if isempty(idx)
  equations(end + 1) = struct('order', order, 'func', {{func}}, 'minus', {{false}});
else
  equations(idx).func{end + 1} = func;
  equations(idx).minus{end + 1} = false;
end

% 更新式方程式の特定次数に負の項を追加
function equations = eqs_minus(equations, order, func)
orders = arrayfun(@(s) s.order, equations);
idx = find(orders == order, 1);
if isempty(idx)
  equations(end + 1) = struct('order', order, 'func', {{func}}, 'minus', {{true}});
else
  equations(idx).func{end + 1} = func;
  equations(idx).minus{end + 1} = true;
end

function eqs_print(equations)
orders = arrayfun(@(s) s.order, equations);
[~, idx] = sort(orders);
for d = idx
  e = equations(d);
  fprintf('  %d order: \n', e.order);
  for n = 1:numel(e.func)
    if e.minus{n}, sig = '-'; else, sig = ''; end
    fprintf('    %s %s\n', sig, func2str(e.func{n}));
  end
end

% 更新式の生成
function [result, str] = eqs_compile(equations, nargs)
orders = arrayfun(@(s) s.order, equations);
min_order = min(orders);
max_order = max(orders);
n = max_order - min_order;
if n == 0
  error('Update function is zeroth order.'); end
c = equations(orders == min_order);
b = equations(orders == min_order + 1);
a = equations(orders == max_order);
if n == 1 && numel(orders) == 2
  % ax + c = 0 -> x = -c/a
  [cf, cfnargs] = sum_func_minus(c.func, c.minus);
  [af, afnargs] = sum_func(a.func, a.minus);
  str = sprintf('(%s)x + %s = 0', func2str(af), func2str(cf));
  result = eval(sprintf('@(%s)cf(%s)./af(%s)', arg_str(nargs), arg_str(cfnargs), arg_str(afnargs)));
elseif isempty(b) && numel(orders) == 2
  % ax^n + c = 0 -> x = (-c/a).^(1/n)
  [cf, cfnargs] = sum_func_minus(c.func, c.minus);
  [af, afnargs] = sum_func(a.func, a.minus);
  str = sprintf('(%s)x^%g + %s = 0', func2str(af), n, func2str(cf));
  result = eval(sprintf('@(%s)(cf(%s)./af(%s)).^(%.16e)', arg_str(nargs), arg_str(cfnargs), arg_str(afnargs), 1 / n));
elseif any(floor(orders) ~= orders)
  % 非整数次を含む上に 3 項以上からなる
  error('The specified equation cannot be solved analytically');
elseif n == 2
  % ax^2 + bx + c = 0
  [af, afnargs] = sum_func(a.func, a.minus);
  [bf, bfnargs] = sum_func(b.func, b.minus);
  [cf, cfnargs] = sum_func(c.func, c.minus);
  str = sprintf('(%s)x^2 + (%s)x + %s = 0', func2str(af), func2str(bf), func2str(cf));
  result = eval(sprintf('@(%s)calceq210(af(%s),bf(%s),cf(%s))',...
    arg_str(nargs), arg_str(afnargs), arg_str(bfnargs), arg_str(cfnargs)));
elseif floor(orders) == orders
  % 整数係数多項式
  % めっちゃ不安定なので推奨されない
  warning('Process will be unstable because the specified equation(s) is high-order / complex polynomial')
  sumfuncs = arrayfun(@(e) sum_func(e.func, e.minus), equations, 'un', 0);
  str = '';
  f = sprintf('@(%s)realmaxroot({', arg_str(nargs));
  for o = n:-1:0
    idx = find(orders == o + min_order, 1);
    if isempty(idx)
      f = sprintf('%s 0', f);
    else
      f = sprintf('%s sumfuncs{%d}(%s)', f, idx, arg_str(nargin(sumfuncs{idx})));
      if o == n
        str = sprintf('(%s)x^%d', func2str(sumfuncs{idx}), o);
      elseif o == 0
        str = sprintf('%s + %s = 0', str, func2str(sumfuncs{idx}));
      else
        str = sprintf('%s + (%s)x^%d', str, func2str(sumfuncs{idx}), n);
      end
    end
  end
  f = sprintf('%s})', f);
  result = eval(f);
else
  error('The specified equation cannot be solved analytically');
end

% 2 次方程式のでかいほうの解
function x = calceq210(a, b, c)
x = (sqrt(b.^2 - 4 * a .* c) - b) ./ a / 2;
% b > 0 のときは精度が落ちるので逆有理化する
flag = b > 0;
x2 = -2 * c ./ (b + sqrt(b.^2 - 4 * a .* c));
x(flag) = x2(flag);


% roots を使って最大の実数解を求める
function x = realmaxroot(p)
sz = max(cell2mat(cellfun(@size, p(:), 'un', 0)), [], 1);
x = arrayfun(@(x, y) realmaxroot_part2(cellfun(@(e) realmaxroot_part1(e, x, y), p)),...
  repmat((1:sz(1))', 1, sz(2)), repmat((1:sz(2)), sz(1), 1));

function s = realmaxroot_part1(s, x, y)
if numel(s) == 1
  return; end
if size(s, 1) == 1
  s = s(1, y);
elseif size(s, 2) == 1
  s = s(x, 1);
else
  s = s(x, y);
end

function x = realmaxroot_part2(p)
x = roots(p);
x = real(x(abs(imag(x)) < eps));
x = max(x);


% 関数を合計する関数を返す
function [func, maxnargs] = sum_func(eqs, minus)
if numel(eqs) == 1 && ~minus{1}
  func = eqs{1};
  maxnargs = nargin(func);
  return;
end

nargs = cellfun(@nargin, eqs);
maxnargs = max(nargs);
str = sprintf('@(%s)', arg_str(maxnargs));
for idx = 1:numel(eqs)
  if minus{idx}
    sig = '-';
  elseif idx == 1
    sig = '';
  else
    sig = '+';
  end
  eval(sprintf('eqs%d = eqs{%d};', idx, idx));
  str = sprintf('%s%seqs%d(%s)', str, sig, idx, arg_str(nargs(idx)));
end
func = eval(str);

% 関数を合計して正負を反転する関数を返す
function [func, maxnargs] = sum_func_minus(eqs, minus)
if numel(eqs) == 1 && minus{1}
  func = eqs{1};
  maxnargs = nargin(func);
  return;
end

nargs = cellfun(@nargin, eqs);
maxnargs = max(nargs);
str = sprintf('@(%s)', arg_str(maxnargs));
for idx = 1:numel(eqs)
  if ~minus{idx}
    sig = '-';
  elseif idx == 1
    sig = '';
  else
    sig = '+';
  end
  eval(sprintf('eqs%d = eqs{%d};', idx, idx));
  str = sprintf('%s%seqs%d(%s)', str, sig, idx, arg_str(nargs(idx)));
end
func = eval(str);

% 引数文字列を生成する
function str = arg_str(arity)
if arity == 0
  str = '';
else
  str = 'arg1';
  for idx = 2:arity
    str = sprintf('%s, arg%d', str, idx); end
end
