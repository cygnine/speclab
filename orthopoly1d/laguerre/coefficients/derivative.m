function[eta] = derivative(n,alpha,varargin)
% [eta] = derivative(n,alpha,{normalization='normal',shift=0, scale=1})
%
%     Computes the connection coefficients between d/dr P^(alpha) and
%     P^(alpha+1). I.e., given the vector n of whole-number indices, and valid
%     Laguerre class parameters alpha, eta is a length(n) vector, where each
%     coefficient is defined by the relation for each q in length(n):
%
%     d/dr*P^(alpha)_n(q) = eta(q)*P^(alpha+1)_(n(q)-1)
%
%     The optional input normalization defines the normalization of the
%     polynomials P, which affects the values of the parameters eta.

persistent defaults
if isempty(defaults)
  from speclab.orthopoly1d.laguerre import defaults
end

opt = defaults(varargin{:});
n = n(:);
N = length(n);

eta = zeros([N,1]);

if strcmpi(opt.normalization,'normal')
  n_is_0 = (n==0);
  n_large = ~(n_is_0);
  temp = n(n_large);
  eta(n_large) = sqrt(temp);
  eta = eta/opt.scale;
else
  fprintf('Error: normalization %s not supported', opt.normalization);
  eta = false;
  return
end
