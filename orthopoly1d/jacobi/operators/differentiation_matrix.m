function[v] = differentiation_matrix(varargin);
% [V] = DIFFERENTIATION_MATRIX({X=GQ(N), N=LENGTH(X), ALPHA=-1/2,BETA=-1/2,SHIFT=0,SCALE=1,NORMALIZATION='NORMAL'});
%
%     Returns the nodal differentiation matrix for the Jacobi polynomials.
%     Either the input X or the input N is required. If e.g. Gauss-Radau is
%     Lobatto points are desired for X, form them externally and feed them in.
%     This function uses direct matrix multiplication to form the differentation
%     matrix. 

global handles;
jac = handles.speclab.orthopoly1d.jacobi;
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

v = handles.common.make_vandermonde(x,n,jac.eval.eval_jacobi_poly,opt);

opt.d = 1;
dv = handles.common.make_vandermonde(x,n,jac.eval.eval_jacobi_poly,opt);

v = dv*inv(v);
