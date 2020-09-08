function writef32(filename, data)
% writef32 - Syntax sugar of writebin(filename, data, 'float32')

srkwii.io.writebin(filename, data, 'float32');
