function writenpy(filename, data, varargin)
% writenpy - saves data into numpy-formatted file
%   writenpy(filename, data)
%   writenpy(___, Name, Value)
%
% Options:
%   Transpose: whether swap 1st and 2nd dimension (default: true)

parser = inputParser;
addRequired(parser, 'filename', @(a) ischar(a) || (isstring(a) && numel(a) == 1));
addRequired(parser, 'data');
addOptional(parser, 'Transpose', true, @isscalar);
parse(parser, input1, input2, varargin{:});

fid = fopen(filename, 'wb');
% magic (\x93NUMPY), version (3.0)
fwrite(fid, [0x93 0x4e 0x55 0x4d 0x50 0x59 0x03 0x00], 'uint8');

[descr, config] = npy_descr(data);
config.transpose = parser.Results.Transpose;

% build header
header = ['{''descr'':''' descr ''',''fortran_order'':False,''shape'':' npy_shape(data, config) '}'];
header_bytes = unicode2native(header, 'UTF-8');
padding = 64 - mod(12 + numel(header_bytes), 64);
fwrite(fid, padding + numel(header_bytes), 'uint32', 0, 'l');
fprintf(fid, '%s%*s', header, padding, newline);

% write data
npy_write(fid, data, config);

fclose(fid);


function [descr, config] = npy_descr(data)
config.numeric_complex = false;

% detect endianness
[~, ~, endian] = computer;
endian_dscr = '<';
if endian == 'B'
  endian_dscr = '>'; end
% detect endianness for string
endian = unicode2native('a', 'UTF-32');
endian_str_dscr = '<';
if endian(1) == 0
  endian_str_dscr = '>'; end

switch class(data)
  case 'double'
    if isreal(data)
      descr = [endian_dscr 'f8'];
    else
      config.numeric_complex = true;
      descr = [endian_dscr 'c16'];
    end
  case 'single'
    if isreal(data)
      descr = [endian_dscr 'f4'];
    else
      config.numeric_complex = true;
      descr = [endian_dscr 'c8'];
    end
  case 'int8'
    descr = 'i1';
  case 'int16'
    descr = [endian_dscr 'i2'];
  case 'int32'
    descr = [endian_dscr 'i4'];
  case 'int64'
    descr = [endian_dscr 'i8'];
  case 'uint8'
    descr = 'u1';
  case 'uint16'
    descr = [endian_dscr 'u2'];
  case 'uint32'
    descr = [endian_dscr 'u4'];
  case 'uint64'
    descr = [endian_dscr 'u8'];
%   case 'datetime'
%   case 'duration'
  case 'string'
    str_maxlen = arrayfun(@strlength, data);
    config.str_maxlen = max(str_maxlen(:));
    descr = [endian_str_dscr 'U' num2str(config.str_maxlen)];
  case 'char'
    config.str_maxlen = size(data, 2);
    descr = [endian_str_dscr 'U' num2str(config.str_maxlen)];
  case 'cell'
    if iscellstr(data) %#ok<ISCLSTR>
      str_maxlen = cellfun(@numel, data);
      config.str_maxlen = max(str_maxlen(:));
      descr = [endian_str_dscr 'U' num2str(config.str_maxlen)];
    else
      error('writenpy not available for cell arrays (except cellstr)');
    end
  otherwise
    error('unsupported type: %s', class(data));
end


function shape = npy_shape(data, config)
sz = size(data);

% sz(2) denots strlen for char array
if ischar(data)
  sz(2) = 1; end

sz = flip(sz);
if config.transpose
  sz(end - 1) = size(data, 1);
  sz(end) = size(data, 2);
end

% omit unnecessary ones
n = find(sz ~= 1, 1);
sz = sz(n:end);

% no elements or one element
if isempty(sz)
  shape = '(1,)';
  return
elseif sum(sz) == 0
  shape = '(0,)';
  return
end

% build string
shape = '(';
for idx = 1:numel(sz)
  if idx > 1
    shape = [shape ',']; end %#ok<AGROW>
  shape = [shape num2str(sz(idx))]; %#ok<AGROW>
end
if numel(sz) == 1
  shape = [shape ',)'];
else
  shape = [shape ')'];
end


function npy_write(fid, data, config)
if config.transpose || ischar(data)
  d = 1:ndims(data);
  d(1:2) = [2 1];
  data = permute(data, d);
end

if isnumeric(data)
  dtype = class(data);
  if config.numeric_complex
    data = [real(data(:))'; imag(data(:))']; end
  fwrite(fid, data, dtype);
  return
end
switch class(data)
  case 'string'
    arrayfun(@(e) npy_writestr(fid, e, config.maxlen), data);
  case 'char'
    r = unicode2native(data(:)', 'UTF-32');
    fwrite(fid, r(5:end));
  case 'cell'
    if iscellstr(data) %#ok<ISCLSTR>
      cellfun(@(e) npy_writestr(fid, e, config.maxlen), data);
    else
      error('writenpy not available for cell arrays (except cellstr)');
    end
end


function npy_writestr(fid, str, len)
r = zeros(1, len * 4, 'uint8');
s = unicode2native(str, 'UTF-32');
if numel(s) > 4
  l = numel(s) - 4;
  r(1:l) = s(5:end);
end
fwrite(fid, r);
