function[eta] = derivative(n,alpha,beta,varargin)
% [ETA] = DERIVATIVE(N,ALPHA,BETA,{NORMALIZATION='normal'})
%
%     Computes the connection coefficients between d/dr P^(ALPHA,BETA) and
%     P^(ALPHA+1,BETA+1). I.e., given the vector N of whole-number indices, and
%     valid Jacobi class parameters ALPHA and BETA, ETA is a length(N)
%     vector, where each coefficient is defined by the relation for each q in
%     length(n):
%
%     d/dr*P^(alpha,beta)_n(q) = eta(q)*P^(alpha+1,beta+1)_(n(q)-1)
%
%     The optional input NORMALIZATION defines the normalization of the
%     polynomials P, which affects the values of the parameters ETA.

global handles;
coeffs = handles.speclab.OrthogonalPolynomial1D.jacobi.coefficients;
opt = handles.common.InputSchema({'normalization'}, {'normal'}, [],varargin{:});
n = n(:);
N = length(n);

eta = zeros([N,1]);

if strcmpi(opt.normalization,'normal')
  n_is_0 = (n==0);
  n_large = ~(n_is_0);
  temp = n(n_large);
  eta(n_large) = sqrt(temp*(temp+alpha+beta+1));
else
  fprintf('Error: normalization %s not supported', opt.normalization);
  mu = false;
  return
end
