function[p] = eval_hermite_poly(x,n,varargin)
% eval_hermite_poly -- evaluate Hermite polynomials
%
% [p] = eval_hermite_poly(x,n,{mu=0,d=0,shift=0,scale=1})
%
%     Evaluates the nth generalized Hermite polynomials at the locations x.
%     This function is vectorized in x and n. The interval of approximation is
%     scaled to interval*scale + shift. The Jacobian resulting from this affine
%     transform is built into the weight function.

persistent opoly hermite
if isempty(opoly)
  from speclab import orthopoly1d as opoly
  from speclab.orthopoly1d import hermite 
end

%global packages;
%opoly = packages.speclab.orthopoly1d;
%hermite = opoly.hermite;
opt = hermite.defaults(varargin{:});

N = max(n)+2;

[a,b] = hermite.coefficients.recurrence(N+1,opt);

if any(strcmpi(opt.normalization,{'normal', 'monic'}))
  p = opoly.eval_polynomial(x,a,b,n,opt);
else
  error('Normalization type not supported');
end

switch opt.weight_normalization
case 'probability'
  p = p*sqrt(b(1));
end
