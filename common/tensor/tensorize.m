function[X] = tensorize(x,dim)
% tensorize -- transforms 1D vector into dim-D array
%
% X = tensorize(x,dim)
%
%     Given a length-N vector x, produces an array X of size (N^dim x dim). The
%     vector x is tensorized over dim dimensions; the resulting N^dim data
%     points in R^dim are ordered (according to Matlab built-in ndgrid ordering)
%     and returned. Being multidimensional, each data point has dim components.

x = x(:);
N = length(x);
Nx = N^dim;

X = zeros([Nx dim]);

xs = cell([dim 1]);
for q = 1:dim
  xs{q} = x;
end

[xs{:}] = ndgrid(xs{:});

for q = 1:dim
  X(:,q) = xs{q}(:);
end
