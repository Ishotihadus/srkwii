function output = cnd(condition, truefunc, falsefunc, varargin)
% cnd - condition(arg1, ..., argn) ? truefunc(arg1, ..., argn) : falsefunc(arg1, ..., argn)
%   output = cnd(condition, truefunc, falsefunc, arg1, ..., argn)

if condition(varargin{:})
  output = truefunc(varargin{:});
else
  output = falsefunc(varargin{:});
end
