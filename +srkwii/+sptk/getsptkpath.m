function path = getsptkpath(command)
% getsptkpath - Get SPTK binary path
%   path = getsptkpath
%   path = getsptkpath(command)
%
% Input:
%   command (string, char[])

global SrkwiiSptkPath

if isempty(SrkwiiSptkPath)
  srkwii.sptk.setsptkpath('/opt/SPTK-3.11-double/bin'); end
path = SrkwiiSptkPath;

if nargin >= 1
  path = fullfile(path, command); end
