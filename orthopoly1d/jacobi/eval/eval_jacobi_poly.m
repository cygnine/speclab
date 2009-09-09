function[p] = eval_jacobi_poly(x,n,varargin)
% [p] = eval_jacobi_poly(x,n,{alpha=-1/2,beta=-1/2,d=0,shift=0,scale=1})
%
%     Evaluates the nth generalized Jacobi polynomials at the locations x. 
%     The weight function is (1-1/scale*(x-shift))^alpha*(1+1/scale*(x-shift))^beta
%     Is vectorized in x and n.

global handles;
opoly = handles.speclab.orthopoly1d;
jac = opoly.jacobi;
opt = jac.defaults(varargin{:});

N = max(n)+2;

[a,b] = jac.coefficients.recurrence(N+1,opt);

if any(strcmpi(opt.normalization,{'normal', 'monic'}))
  p = opoly.eval_polynomial(x,a,b,n,opt);
else
  error('Normalization type not supported');
end
