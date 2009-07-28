function[p] = eval_jacobi_poly(x,n,varargin)
% [P] = EVAL_JACOBI_POLY(X,N,{ALPHA=-1/2,BETA=-1/2,SHIFT=0,SCALE=1})
%
%     Evaluates the nth generalized Jacobi polynomials at the locations x. 
%     The weight function is (1-1/scale*(x-shift))^alpha*(1+1/scale*(x-shift))^beta
%     Is vectorized in x and n.

global handles;
opoly = handles.speclab.OrthogonalPolynomial1D;
jac = opoly.jacobi;
opt = jac.defaults(varargin{:});

N = max(n)+2;

[a,b] = jac.recurrence(N+1,opt);

if strcmpi(opt.normalization,'normal')
  p = opoly.eval_polynomial(x,a,b,n,opt);
end
