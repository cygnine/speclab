function[X] = tensorize_grid1d(x,dim)
% tensorize_grid1d -- transforms a 1D vector into a dim-D array
%
% X = tensorize_grid1d(x,dim)
%
%     Given a length-N vector x, produces an array X of size (N^dim x dim). The
%     vector x is tensorized over dim dimensions; the resulting N^dim data
%     points in R^dim are ordered (according to Matlab built-in ndgrid ordering)
%     and returned. Being multidimensional, each data point has dim components.

persistent tensorize
if isempty(tensorize)
  from speclab.common.tensor import tensorize
end

X = tensorize(x,dim);
