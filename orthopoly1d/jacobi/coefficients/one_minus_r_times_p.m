function[mu] = one_minus_r_times_p(n,alpha,beta,varargin)
% [mu] = one_minus_r_times_p(n,alpha,beta,{normalization='normal'})
%
%     Computes the connection coefficients between (1-r) x P and P. I.e., given
%     the vector n of whole-number indices, and valid Jacobi class parameters
%     alpha and beta, mu is a length(n) x 2 matrix, where each row contains the
%     coefficients defined by the relation:
%
%     (1-r)*P^(alpha,beta)_n = mu(1)*P^(alpha-1,beta)_n +
%                              mu(2)*P^(alpha-1,beta)_{n+1}
%
%     The optional input normalization defines the normalization of the
%     polynomials P, which affects the values of the parameters mu.

persistent input_schema
if isempty(input_schema)
  from labtools import input_schema
end

opt = input_schema({'normalization','scale'}, {'normal',1}, [],varargin{:});
n = n(:);
N = length(n);

mu = zeros([N,2]);

if strcmpi(opt.normalization,'normal')
  tol = 1e-12;
  if abs(alpha+beta)<tol
    n_is_0 = n==0;
    n_not_0 = ~(n==0);
    mu(n_is_0,1) = sqrt(2*alpha/(alpha+beta+1));
    ns = n(n_not_0);
    mu(n_not_0,1) = sqrt(2*(ns + alpha).*(ns+alpha+beta)./...
                         ((2*ns+alpha+beta).*(2*ns+alpha+beta+1)));
  else
    mu(:,1) = sqrt(2*(n+alpha).*(n+alpha+beta)./...
            ((2*n+alpha+beta).*(2*n+alpha+beta+1)));
  end

  mu(:,2) = -sqrt(2*(n+1).*(n+beta+1)./...
            ((2*n+alpha+beta+1).*(2*n+alpha+beta+2)));

  mu = mu/opt.scale;
else
  fprintf('Error: normalization %s not supported', opt.normalization);
  mu = false;
  return
end
