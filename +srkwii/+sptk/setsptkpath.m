function setsptkpath(path)
% setsptkpath - Set SPTK binary path
%
% Input:
%   path (string, char[])

global SrkwiiSptkPath SrkwiiSptkIsDouble

SrkwiiSptkPath = path;

outfile = tempname;
builder = srkwii.sptk.CommandBuilder('step');
builder.AddOption('l', 'int', 1);
builder.SetOutput(outfile);
builder.Run;
output = srkwii.io.readu32(outfile);
SrkwiiSptkIsDouble = numel(output) == 2;
delete(outfile);
