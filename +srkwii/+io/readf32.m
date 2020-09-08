function output = readf32(filename, varargin)
% readf32 - Syntax sugar of readbin(filename, 'float32', dim1, ..., dimn)

output = srkwii.io.readbin(filename, 'float32', varargin{:});
