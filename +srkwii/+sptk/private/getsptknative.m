function out = getsptknative
% getsptknative - Get SPTK native type

global SrkwiiSptkIsDouble
out = 'float32';
if SrkwiiSptkIsDouble
  out = 'float64'; end
