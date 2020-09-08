function [index, str] = selectstr(query, list, varargin)
% selectstr - validate and select str
%   index = selectstr(query, list)
%   index = selectstr(query, list, argIdx)
%   index = selectstr(query, list, funcName)
%   index = selectstr(query, list, funcName, varName)
%   index = selectstr(query, list, funcName, varName, argIdx)
%   [index, str] = selectstr(___)

if isnumeric(query)
  index = query;
  if nargout > 1
    str = list{index + 1};
    if iscell(str)
      str = str{1}; end
  end
  return
end

list_expand = cellfun(@flatcell, list, 'un', 0);
list_expand_idx = cell2mat(...
  arrayfun(@(e, n) repmat(n, 1, e), cellfun(@numel, list_expand(:)'), 1:numel(list_expand), 'un', 0));
list_expand = [list_expand{:}];
str = validatestring(query, list_expand, varargin{:});
matches = cellfun(@(e) strcmp(e, str), list_expand);
index = list_expand_idx(matches) - 1;


function c = flatcell(c)
if iscell(c)
  c = c(:)';
else
  c = {c};
end
