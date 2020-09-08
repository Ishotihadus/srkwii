classdef HandleStruct < matlab.mixin.Copyable
  properties (SetAccess = private)
    Data
  end

  methods
    function obj = HandleStruct(varargin)
      obj.Data = struct(varargin{:});
    end

    function varargout = subsref(obj, s)
      varargout = cell(1, max(1, nargout));
      [varargout{:}] = builtin('subsref', obj.Data, s);
    end
    function obj = subsasgn(obj, s, a)
      obj.Data = builtin('subsasgn', obj.Data, s, a); end

    function n = numArgumentsFromSubscript(obj, s, indexingContext)
      n = numArgumentsFromSubscript(obj.Data, s, indexingContext); end
    function sz = size(obj, varargin)
      sz = size(obj.Data, varargin{:}); end
    function sz = numel(obj)
      sz = numel(obj.Data); end

    function fields = fieldnames(obj, varargin)
      fields = fieldnames(obj.Data, varargin{:}); end
    function value = getfield(obj, varargin)
      value = getfield(obj.Data, varargin{:}); end
    function TF = isfield(obj, varargin)
      TF = builtin('isfield', obj.Data, varargin{:}); end
    function T = isstruct(obj)
      T = true; end
    function obj = orderfields(obj, varargin)
      obj.Data = orderfields(obj.Data, varargin{:}); end
    function obj = rmfield(obj, varargin)
      obj.Data = rmfield(obj.Data, varargin{:}); end
    function obj = setfield(obj, varargin)
      obj.Data = setfield(obj.Data, varargin{:}); end

    function disp(obj)
      fprintf('  srkwii.util.HandleStruct:\n\n');
      disp(obj.Data);
    end
  end
end
