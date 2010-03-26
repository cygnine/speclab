function[coeffs] = exp_gegenbauer_expansion(theta,n,varargin)
% exp_gegenbauer_expansion -- Gegenbauer polynomial expansion cofficients
%
% coeffs = exp_gegenbauer_expansion(theta,n,{lambda=0})
%
%     Computes projection of the function exp(theta * x) over x \in [-1,1] onto
%     the n-th Gegenbauer polynomial in x of class lambda. This is given by an
%     analytic formula involving Bessel functions. The output coeffs is a
%     length(theta) x length(n) matrix, with the coefficients. The Gegenbauer
%     polynomials are assumed to be orthonormal in the native L^2 weight.
%
%     Lambda can take any value greater than -1/2; the code is not vectorized in
%     lambda.

persistent strict_inputs spdiag recurrence
if isempty(strict_inputs)
  from labtools import strict_inputs spdiag
  from speclab.orthopoly.jacobi.coefficients import recurrence
end

opt = strict_inputs({'lambda'}, {0}, [], varargin{:});

theta = theta(:);
n = n(:);

% For theta==0, we'll deal with it at the end...right now just make a note of it
% and set it to some nonzero value.
theta_zero = theta==0;
theta(theta_zero) = 1;

assert(all(n>=0), 'Error: the coefficient indices n must be non-negative');
assert(length(opt.lambda)==1, 'Error: lambda must be a scalar');
assert(opt.lambda>-1/2, 'Error: lambda must be greater than -1/2');

sign_theta = sign(theta);
theta = abs(theta);

coeffs = besseli(n.'+opt.lambda+1/2, theta);

n_factors = sqrt(2*pi)*sqrt(n+opt.lambda+1/2);
n_factors = n_factors.*sqrt(exp(gammaln(n+2*opt.lambda+1) - gammaln(n+1)));

theta_factors = 1./theta.^(opt.lambda+1/2);

coeffs = spdiag(theta_factors)*coeffs*spdiag(n_factors);

signs = repmat(sign_theta, [1 length(n)]);
signs = signs.^repmat(n.', [length(theta) 1]);

coeffs = coeffs.*signs;

% If theta==0, then the coefficient n=0 is the only nonzero one
if any(theta_zero);
  coeffs(theta_zero,:) = 0;
  [a,b] = recurrence(0, 'alpha', opt.lambda, 'beta', opt.lambda);
  coeffs(theta_zero, n==0) = 1/sqrt(b);
end
