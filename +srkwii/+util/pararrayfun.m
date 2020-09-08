function varargout = pararrayfun(func, varargin)
% pararrayfun - arrayfun with parfor
%   B = pararrayfun(func, A1, ..., An)
%   B = pararrayfun(___, Name, Value)
%   [B1, ..., Bm] = pararrayfun(___)

narginchk(2, inf);

for idx = 1:numel(varargin)
  if ischar(varargin{idx}) || isstring(varargin{idx})
    break; end
  varargin{idx} = num2cell(varargin{idx});
end

varargout = cell(1, nargout);
[varargout{:}] = srkwii.util.parcellfun(func, varargin{:});
