function[x,w] = gauss_quadrature(a,b,N,varargin);
% [x,w] = gauss_quadrature(a,b,n,{dim=1});
%
%     Returns the N-point Gaussian polynomial quadrature for the orthogonal
%     polynomials corresponding to the recurrence coefficients A and B.
%
%     Returns a tensorized grid with N degrees of freedom in each dimension.
%     I.e. returns an Nx x dim matrix, where Nx = N^dim.

persistent strict_inputs tensorize
if isempty(strict_inputs)
  from labtools import strict_inputs
  from speclab.common.tensor import tensorize
end

opt = strict_inputs({'dim'}, {1}, [], varargin{:});

a = a(1:N);
b = b(1:N);

a = a(:);
b = b(:);

J = spdiags([sqrt([b(2:end);0]) a sqrt(b)], -1:1, N, N);

[v,d] = eig(full(J));
x = diag(d);
w = b(1)*v(1,:).^2.';

if opt.dim>1
  x = tensorize(x, opt.dim);
  w = tensorize(w, opt.dim);
  w = prod(w, 2);
end
