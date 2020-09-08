function varargout = parcellfun(func, varargin)
% parcellfun - cellfun with parfor
%   A = parcellfun(func, C1, ..., Cn)
%   A = parcellfun(___, Name, Value)
%   [A1, ..., Am] = parcellfun(___)

narginchk(2, inf);

for numcells = 1:numel(varargin)
  if ischar(varargin{numcells}) || isstring(varargin{numcells})
    numcells = numcells - 1;
    break
  elseif iscell(varargin{numcells})
    if numcells > 1 && ~all(size(varargin{numcells}) == size(varargin{1}))
      error('All of the input arguments must be of the same size and shape.'); end
  else
    error('Input #%d expected to be a cell array, was %s instead.', numcells + 1, class(varargin{numcells}))
  end
end

parser = inputParser;
addOptional(parser, 'UniformOutput', true, @(a) isscalar(a));
parse(parser, varargin{numcells + 1:end});

shape = size(varargin{1});
count = numel(varargin{1});

argslist = cellfun(...
  @(idx) cellfun(@(a) a{idx}, varargin(1:numcells), 'UniformOutput', false), num2cell(1:count), 'UniformOutput', false);

numoutput = nargout;
if nargout == 0
  func_nargout = nargout(func);
  if func_nargout < -1 || func_nargout > 0
    numoutput = 1; end
end

if numoutput == 0
  parfor idx = 1:count
    args = argslist{idx};
    func(args{:});
  end
  return
elseif numoutput == 1
  output = cell(count, 1);
  parfor idx = 1:count
    args = argslist{idx};
    output{idx} = func(args{:});
  end
  output = reshape(output, shape);
  varargout{1} = output;
else
  output = cell(count, numoutput);
  parfor idx = 1:count
    args = argslist{idx};
    tmpoutput = cell(1, numoutput);
    [tmpoutput{:}] = func(args{:});
    output(idx, :) = tmpoutput;
  end
  output = cellfun(@(idx) reshape(output(:, idx), shape), num2cell(1:numoutput), 'UniformOutput', false);
  varargout = output;
end

if parser.Results.UniformOutput
  cellfun(@(idx) checkun(varargout{idx}, idx), num2cell(1:numel(varargin)));
  varargout = cellfun(@(a) cell2mat(a), varargout, 'UniformOutput', false);
end


function checkun(input, id)
k = find(~cellfun(@isscalar, input), 1);
if ~isempty(k)
  error('Non-scalar in Uniform output, at index %d, output %d.\nSet ''UniformOutput'' to false.', k, id); end
