function[X] = subarray_tensorize(x, dim, varargin)
% subarray_tensorize -- Returns a subarray of a tensorized vector
%
% X = subarray_tensorize(x, dim, {index=[]})
%
%     Returns a subarray of tensorize_grid1d. In high-dimensional settings, it
%     is not possible to generate an array of size (length(x)^dim) x dim, which
%     would be a dim-dimensional tensor product of x. This function returns a
%     submatrix of the above. In particular, if Y = tensorize_grid1d(x,dim),
%     then X = Y(index,:). If `index' is empty, the entire array is returned. To
%     generate the linear indices from dim-dimensional indexing, see Matlab's
%     sub2ind and ind2sub.
%
%     This function returns the subarray without internally generating the full
%     array, thus preventing memory issues if `index' is a sufficiently small
%     subset. 

persistent tensorize strict_inputs
if isempty(tensorize)
  from speclab.grids.tensor import tensorize_grid1d as tensorize
  from labtools import strict_inputs
end

opt = strict_inputs({'index'}, {[]}, [], varargin{:});

if isempty(opt.index)
  X = tensorize(x,dim);
  return
end

X = zeros([length(opt.index) dim]);

Y = cell([dim 1]);
[Y{:}] = ind2sub(length(x)*ones([1 dim]), opt.index);

for q = 1:dim
  X(:,q) = x(Y{q});
end
