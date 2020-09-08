function list = readlist(filename)
% readlist - Read text and store into a vertical cell
%
% Output:
%   list [cellstr{n, 1}]

fileID = fopen(filename);
list = textscan(fileID, '%s');
list = list{1}(:);
fclose(fileID);
