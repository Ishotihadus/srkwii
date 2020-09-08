function varargout = struct2cell(S, varargin)
% struct2cell - convert structure to cell array
%   [C] = struct2cell(S)
%   [C1, ..., CN] = struct2cell(S, field1, ..., fieldN)

narginchk(1, inf);

switch nargin
  case 1
    varargout{1} = struct2cell(S);
  case 2
    varargout{1} = arrayfun(@(s) s.(varargin{1}), S,  'un', 0);
  otherwise
    varargout = cellfun(@(f) arrayfun(@(s) s.(f), S,  'un', 0), varargin,  'un', 0);
end
