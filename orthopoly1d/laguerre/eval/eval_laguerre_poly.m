function[p] = eval_laguerre_poly(x,n,varargin)
% eval_laguerre_poly -- evaluate Laguerre polynomials
%
% [p] = eval_laguerre_poly(x,n,{alpha=0,d=0,shift=0,scale=1})
%
%     Evaluates the nth generalized Laguerre polynomials at the locations x.
%     This function is vectorized in x and n. The interval of approximation is
%     scaled to interval*scale + shift. The Jacobian resulting from this affine
%     transform is built into the weight function.

global packages;
opoly = packages.speclab.orthopoly1d;
laguerre = opoly.laguerre;
opt = laguerre.defaults(varargin{:});

N = max(n)+2;

[a,b] = laguerre.coefficients.recurrence(N+1,opt);

if any(strcmpi(opt.normalization,{'normal', 'monic'}))
  p = opoly.eval_polynomial(x,a,b,n,opt);
else
  error('Normalization type not supported');
end
