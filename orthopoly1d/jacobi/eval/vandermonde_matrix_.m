function[v] = vandermonde_matrix(varargin);
% [v] = vandermonde_matrix({x=gq(n), n=length(x), alpha=-1/2,beta=-1/2,shift=0,scale=1,normalization='normal'});
%
%     Returns the Vandermonde matrix for the Jacobi polynomials. Either the
%     input x or the input n is required. If e.g. Gauss-Radau is Lobatto
%     points are desired for x, form them externally and feed them in.

persistent defaults make_vandermonde gq
if isempty(defaults)
  from speclab.orthopoly1d.jacobi import defaults
  from speclab.orthopoly1d.jacobi.quad import gauss_quadrature as gq
  from labtools import make_vandermonde
end

opt = defaults(varargin{:});

if x ~= false
  x = x(:);
  n = 0:(length(x)-1);

elseif n ~= false
  [x,w] = gq(n,opt);
  n = 0:(n-1);
else
  error('You must input either a nodal vector x or a size n');
end

v = make_vandermonde(x,n,jac.eval.eval_jacobi_poly,opt);
