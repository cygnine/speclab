function[eta] = derivative(n,alpha,beta,varargin)
% [eta] = derivative(n,alpha,beta,{normalization='normal',shift=0, scale=1})
%
%     Computes the connection coefficients between d/dr P^(alpha,beta) and
%     P^(alpha+1,beta+1). I.e., given the vector n of whole-number indices, and
%     valid Jacobi class parameters alpha and beta, eta is a length(n)
%     vector, where each coefficient is defined by the relation for each q in
%     length(n):
%
%     d/dr*P^(alpha,beta)_n(q) = eta(q)*P^(alpha+1,beta+1)_(n(q)-1)
%
%     The optional input normalization defines the normalization of the
%     polynomials P, which affects the values of the parameters eta.

global packages;
coeffs = packages.speclab.orthopoly1d.jacobi.coefficients;
opt = packages.common.input_schema({'normalization','scale'}, {'normal',1}, [],varargin{:});
n = n(:);
N = length(n);

eta = zeros([N,1]);

if strcmpi(opt.normalization,'normal')
  n_is_0 = (n==0);
  n_large = ~(n_is_0);
  temp = n(n_large);
  eta(n_large) = sqrt(temp.*(temp+alpha+beta+1));
  eta = eta/opt.scale;
else
  fprintf('Error: normalization %s not supported', opt.normalization);
  eta = false;
  return
end
