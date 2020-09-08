function output = readu32(filename, varargin)
% readu32 - Syntax sugar of readbin(filename, 'uint32', dim1, ..., dimn)

output = srkwii.io.readbin(filename, 'uint32', varargin{:});
