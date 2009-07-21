function[mu] = one_minus_r_times_p(n,alpha,beta,varargin)
% [MU] = ONE_MINUS_R_TIMES_P(N,ALPHA,BETA,{NORMALIZATION='normal'})
%
%     Computes the connection coefficients between (1-r) x P and P. I.e., given
%     the vector N of whole-number indices, and valid Jacobi class parameters
%     ALPHA and BETA, MU is a length(N) x 2 matrix, where each row contains the
%     coefficients defined by the relation:
%
%     (1-r)*P^(alpha,beta)_n = mu(1)*P^(alpha-1,beta)_n +
%                              mu(2)*P^(alpha-1,beta)_{n+1}
%
%     The optional input NORMALIZATION defines the normalization of the
%     polynomials P, which affects the values of the parameters MU.

global handles;
opt = handles.common.InputSchema({'normalization'}, {'normal'}, [],varargin{:});
n = n(:);
N = length(n);

mu = zeros([N,2]);

switch opt.normalization
case 'normal'
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
otherwise
  fprintf('Error: normalization %s not supported', opt.normalization);
  mu = false;
  return
end
