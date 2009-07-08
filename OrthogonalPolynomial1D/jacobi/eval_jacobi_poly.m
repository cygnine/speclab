function[p] = eval_jacobi_poly(x,n,varargin)
% [P] = EVAL_JACOBI_POLY(X,N,{ALPHA=-1/2,BETA=-1/2,SHIFT=0,SCALE=1})
%
%     Evaluates the nth generalized Jacobi polynomials at the locations x. 
%     The weight function is (1-1/scale*(x-shift))^alpha*(1+1/scale*(x-shift))^beta
%     Is vectorized in x and n.
%
% 20080524: acn

global handles;
opoly = handles.speclab.OrthogonalPolynomial1D;
jac = opoly.jacobi;
sss = handles.speclab.common.standard_scaleshift_1d;
opt = jac.defaults(varargin{:});

N = max(n)+2;

[a,b] = jac.recurrence(N+1,opt);

switch opt.normalization
case 'normal'
  p = opoly.eval_polynomial(x,a,b,n,opt);
end
