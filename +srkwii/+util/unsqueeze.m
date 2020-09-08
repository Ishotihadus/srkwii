function A = unsqueeze(A, dim)
% unsqueeze - returns a new tensor with a dimension of size one inserted at the specified position
%   B = unsqueeze(A, dim)
%
% Input:
%   dim: unsqueezed dimension (default: 1)

if nargin <= 1
  dim = 1; end

if ndims(matrix) >= dim
  sz = size(matrix);
  sz = [sz(1:dim - 1) 1 sz(dim:end)];
  A = reshape(A, sz);
end
