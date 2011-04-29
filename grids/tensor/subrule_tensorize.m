function[X,W] = subrule_tensorize(x, w, dim, varargin)
% subrule_tensorize -- Returns a subrule of a tensorized quadrature rule
%
% X = subrule_tensorize(x, w, dim, {index=[]})
%
%     Returns a subset of a quadrature rule that is a dim-dimensional tensorized
%     copy of the 1D rule (x,w). The motivation for this function is the same as
%     that for subarray_tensorize, and this function merely constructs a
%     quadrature rule subset instead of just the nodes.

persistent subarray_tensorize
if isempty(subarray_tensorize)
  from speclab.grids.tensor import subarray_tensorize
end

X = subarray_tensorize(x,dim, varargin{:});
W = subarray_tensorize(w,dim, varargin{:});
W = prod(W,2);
