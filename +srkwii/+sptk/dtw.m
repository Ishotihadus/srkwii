function [vit, aligned1, aligned2, score] = dtw(input1, input2, varargin)
% dtw - dynamic time warping
%   [vit, aligned1, aligned2, score] = dtw(input1, input2)
%   [vit, aligned1, aligned2, score] = dtw(___, Name, Value)
%
% Options:
%   Norm     : type of norm used for calculation of local cost (1 or 2, default: 2)
%   LocalPath: local path constraint (1–7, default: 5)
%     1:       2:       3:
%       o > o    o > o    o > o
%           ^      / ^      /
%           o    o   o    o
%
%     4:       5:           6:           7:
%       o > o        o > o            o        o > o
%         //       /   / |        / //       //  //
%       o /      o   o   o    o   o /      o / o /
%        /             /           /        /   /
%       o            o            o        o   o

parser = inputParser;
parser.addRequired('input1', @isreal);
parser.addRequired('input2', @isreal);
parser.addOptional('Norm', 2, @(a) isscalar(a) && (a == 1 || a == 2));
parser.addOptional('LocalPath', 5, @(a) isscalar(a) && 1 <= a && a <= 7);
parser.parse(input1, input2, varargin{:});

if size(input1, 1) ~= size(input2, 1)
  error('the length of vectors is not same'); end

dim = size(input1, 1);

infile1 = tempname;
infile2 = tempname;
scorefile = tempname;
vitfile = tempname;

srkwii.io.writebin(infile1, input1, getsptknative);
srkwii.io.writebin(infile2, input2, getsptknative);

builder = srkwii.sptk.CommandBuilder('dtw');
builder.AddOption('l', 'int', dim);
builder.AddOption('n', 'int', parser.Results.Norm);
builder.AddOption('p', 'int', parser.Results.LocalPath);
if nargout >= 4
  builder.AddOption('s', 'string', scorefile); end
builder.AddOption('v', 'string', vitfile);
builder.AddArgument(infile2);
builder.AddArgument(infile1);
builder.Run;

vit = srkwii.io.readbin(vitfile, 'int', 2) + 1;
if nargout >= 2
  aligned1 = input1(:, vit(1, :));
  aligned2 = input2(:, vit(2, :));
end
if nargout >= 4
  score = srkwii.io.readbin(scorefile, getsptknative, 1); end

delete(infile1);
delete(infile2);
if nargout >= 4
  delete(scorefile); end
delete(vitfile);
