function[mu] = one_plus_r_times_p(n,alpha,beta,varargin)
% [mu] = one_plus_r_times_p(n,alpha,beta,{normalization='normal'})
%
%     Computes the connection coefficients between (1+r) x P and P. I.e., given
%     the vector N of whole-number indices, and valid Jacobi class parameters
%     alpha and beta, mu is a length(N) x 2 matrix, where each row contains the
%     coefficients defined by the relation:
%
%     (1+r)*P^(alpha,beta)_n = mu(1)*P^(alpha,beta-1)_n +
%                              mu(2)*P^(alpha,beta-1)_{n+1}
%
%     The optional input normalization defines the normalization of the
%     polynomials P, which affects the values of the parameters mu.

global packages;
coeffs = packages.speclab.orthopoly1d.jacobi.coefficients;
opt = packages.common.input_schema({'normalization','scale'}, {'normal',1}, [],varargin{:});
n = n(:);
N = length(n);

mu = zeros([N,2]);

if strcmpi(opt.normalization,'normal')
  mu = coeffs.one_minus_r_times_p(n,beta,alpha,opt);
  mu(:,2) = -mu(:,2);
else
  fprintf('Error: normalization %s not supported', opt.normalization);
  mu = false;
  return
end
