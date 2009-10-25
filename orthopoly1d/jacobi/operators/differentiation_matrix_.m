function[v] = differentiation_matrix(varargin);
% [v] = differentiation_matrix({x=gq(n), n=length(x), alpha=-1/2,beta=-1/2,shift=0,scale=1,normalization='normal'});
%
%     Returns the nodal differentiation matrix for the Jacobi polynomials.
%     Either the input x or the input n is required. If e.g. Gauss-Radau is
%     Lobatto points are desired for x, form them externally and feed them in.
%     This function uses direct matrix multiplication to form the differentation
%     matrix. 

global packages;
jac = packages.speclab.orthopoly1d.jacobi;
opt = jac.defaults(varargin{:});

if x ~= false
  x = x(:);
  n = 0:(length(x)-1);

elseif n ~= false
  [x,w] = jac.quad.gauss_quadrature(n,opt);
  n = 0:(n-1);
else
  error('You must input either a nodal vector X or a size N');
end

v = packages.labtools.make_vandermonde(x,n,jac.eval.eval_jacobi_poly,opt);

opt.d = 1;
dv = packages.labtools.make_vandermonde(x,n,jac.eval.eval_jacobi_poly,opt);

v = dv*inv(v);
