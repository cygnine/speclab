function[mu] = one_minus_r_squared_times_p(n,alpha,beta,varargin)
% [MU] = ONE_MINUS_R_SQUARED_TIMES_P(N,ALPHA,BETA,{NORMALIZATION='normal',SCALE=1})
%
%     Computes the connection coefficients between (1-r^2) x P and P. I.e.,
%     given the vector N of whole-number indices, and valid Jacobi class
%     parameters ALPHA and BETA, MU is a length(N) x 3 matrix, where each row
%     contains the coefficients defined by the relation:
%
%     (1-r^2)*P^(alpha,beta)_n = mu(1)*P^(alpha-1,beta-1)_n +
%                                mu(2)*P^(alpha-1,beta-1)_{n+1}
%                                mu(3)*P^(alpha-1,beta-1)_{n+2}
%
%     The optional input NORMALIZATION defines the normalization of the
%     polynomials P, which affects the values of the parameters MU.

global handles;
coeffs = handles.speclab.OrthogonalPolynomial1D.jacobi.coefficients;
opt = handles.common.InputSchema({'normalization','scale'}, {'normal',1}, [],varargin{:});
n = n(:);
N = length(n);

mu = zeros([N,3]);

if strcmpi(opt.normalization,'normal')
  n_is_0 = (n==0);
  n_is_1 = (n==1);
  n_large = ~(n_is_0 & n_is_1);
  temp = 2*n(n_large)+alpha+beta;
  mu(n_large,1) = sqrt(4*(n(n_large)+alpha).*(n(n_large)+beta).*(n(n_large)+alpha+beta-1).*(n(n_large)+alpha+beta)./...
                 ((temp-1).*(temp.^2).*(temp+1)));
  mu(n_is_0,1) = sqrt(4*alpha*beta./...
                      ((alpha+beta)*(alpha+beta+1)));
  mu(n_is_1,1) = 2/(alpha+beta+2)*sqrt((alpha+1)*(beta+1)*(alpha+beta)/(alpha+beta+3));
%    epsn[:,0] = sqrt(4*(n+a)*(n+b)*(n+a+b-1)*(n+a+b)/ \
%                        ((2*n+a+b-1)*(2*n+a+b)**2*(2*n+a+b+1)))

  temp = 2*n+alpha+beta;
  mu(:,2) = 2*(alpha-beta)*sqrt((n+1).*(n+alpha+beta))./...
                                (temp.*(temp+2));
%    epsn[:,1] = 2*(alpha-beta)*sqrt((n+1)*(n+a+b))/ \
%                ((2*n+a+b)*(2*n+a+b+2))

  mu(:,3) = -2*sqrt((n+1).*(n+2).*(n+alpha+1).*(n+beta+1)./...
                    ((temp+1).*(temp+2).^2.*(temp+3)));
%    epsn[:,2] = -sqrt(4*(n+1)*(n+2)*(n+a+1)*(n+b+1)/ \
%                           ((2*n+a+b+1)*(2*n+a+b+2)**2*(2*n+a+b+3)))
  mu = mu/opt.scale;
else
  fprintf('Error: normalization %s not supported', opt.normalization);
  mu = false;
  return
end
