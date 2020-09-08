function path = checkcommandpath(path)
% checkcommandpath - Check whether the binary path is valid

if isfolder(path)
  error('`%s`: is a directory', path); end

if ~isfile(path)
  error('`%s`: no such file or directory', path); end

[attrib_status, attrib_values] = fileattrib(path);

if ~attrib_status
  error('`%s`: failed to read attributes (permission denied?)', path); end

if ~attrib_values.UserExecute
  error('`%s`: not executable', path); end
