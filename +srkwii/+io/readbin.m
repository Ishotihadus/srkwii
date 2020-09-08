function output = readbin(filename, precision, varargin)
% readbin - Read SPTK-formatted binary and reshape
%   output = readbin(filename, precision)
%   output = readbin(filename, precision, dim1, ..., dimn)

fileID = fopen(filename);
output = fread(fileID, Inf, precision);
fclose(fileID);

% 空っぽの場合はそのまま返す
if isempty(varargin)
  return; end

% Inf の要素を含む場合は Inf の要素を [] にする
% Inf の要素を含まない場合には最後の要素に [] を追加する
valflag = cellfun(@(s) numel(s) == 1 && s ~= Inf, varargin);
if all(valflag)
  varargin{1, end + 1} = [];
else
  varargin{1, ~valflag} = [];
end

output = reshape(output, varargin{:}, 1);
