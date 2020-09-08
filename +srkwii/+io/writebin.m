function writebin(filename, data, precision)
% writebin - Write SPTK-formatted binary

fileID = fopen(filename, 'w');
fwrite(fileID, data, precision);
fclose(fileID);
