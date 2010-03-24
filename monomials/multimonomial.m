function[fz] = multimonomial(x, n, varargin)
% multinomial -- Evaluates multidimensional monomials
%
% fz = multimonomial(x, n, [dim=1])
%
%     Evaluates multidimensional monomials at the locations x. The input
%     indexing n should be linear indexing; this is then translated into
%     multidimensional indexing using
%     speclab.common.tensor.linear_to_array_indexing.
%
%     The output fz is a length(x(:)) x length(n(:)) array. The optional input
%     "dim" specifies the dimension. In cases when dim>1, then x should have dim
%     columns, and the output fz has size(x,1) rows.
%
%     Example:
%
%       fz = multimonomial([3, 5; 
%                           4, 6;
%                           8, 9;], [0, 1, 2], 'dim', 2);
%
%       returns fz as a 3 x 3 matrix. Column one contains the evaluations of
%       function 0 (constant 1) at the 2D points (3,5), (4, 6), and (8,9).
%       Columns two and three then contain the evaluations of functions 1 and 2
%       (x and y, respectively) at the same locations.

persistent indexing strict_inputs
if isempty(strict_inputs)
  from labtools import strict_inputs
  from speclab.common.tensor import linear_to_array_indexing as indexing
end

opt = strict_inputs({'dim'}, {1}, [], varargin{:});
dim = opt.dim;

indices = indexing(n(:), 'dim', dim);

if dim>1
  assert(size(x,2)==dim, 'Error: x must have the same number of columns as the dimension');

  fz = ones(size(x,1), size(indices,1));
else
  x = x(:);
  fz = ones(length(x), size(indices,1));
end

[M,N] = size(fz);

% With potentially random ordering for indices, I have no idea how to do a
% nested Horner-type algorithm
for q = 1:opt.dim;
  fz = fz.*repmat(x(:,q), [1 N]).^repmat(indices(:,q).', [M 1]);
end
