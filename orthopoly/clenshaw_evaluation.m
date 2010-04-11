function[f] = clenshaw_evaluation(x, coeffs, alpha, beta, options)
% clenshaw_evaluation -- Evaluates an orthogonal polynomial sum
%
% f = clenshaw_evaluation(x, coeffs, alpha, beta, options)
%
%     Evaluates the truncated polynomial sum
%
%       f(x) = \sum_n coeffs(n)*p_{n-1}(x),
%
%     where the polynomials p_n are generated from a three-term recurrence given
%     by the alpha and beta coefficients. A test is performed to ensure that
%     enough coefficients are given to match the length of coeffs. This function
%     is not vectorized if x and coeffs are matrices, but it can surely be done
%     (TODO).
%
%     'options' is a struct with the following fields:
%
%           options.normalization ({'normal'}, 'monic') -- The normalization
%           used for the polynomials p_n. See speclab.orthopoly.eval_polynomial.
%
%           options.shift ({0}) -- The affine shift p(x) <---- p(x-shift)
%
%           options.scale ({1}) -- The affine scale p(x) <---- p(x/scale);

persistent sss
if isempty(sss)
  from speclab.common import standard_scaleshift_1d as sss
end

x = sss(x, options);

x = x(:);
N = length(coeffs);
Nx = length(x);
A = length(alpha);
B = length(beta);

switch options.normalization
case 'monic'
  if (A<N) | (B<N)
    error('For a monic formulation, I need at least N recurrence coefficients');
  end

  storage = zeros([Nx 1]);
  bb = zeros([Nx 1]);  % b_N
  b = coeffs(N)*ones([Nx 1]); % b_{N-1}

  for q = (N-1):-1:2
    storage = bb;
    bb = b;
    b = coeffs(q) - beta(q+1)*storage + (x-alpha(q)).*b;
  end

  % f = b1*p1 + a0 - b2*beta1
  f = b.*((x-alpha(1))) + coeffs(1) - bb*beta(2);
case 'normal'
  if (A<N+1) | (B<N+1)
    error('For an orthonormal formulation I need at least N+1 recurrence coefficients');
  end
  beta = sqrt(beta);

  storage = zeros([Nx 1]);
  bb = zeros([Nx 1]);
  b = coeffs(N)*ones([Nx 1]);

  for q = (N-1):-1:2
    storage = bb;
    bb = b;
    b = coeffs(q) - beta(q+1)/beta(q+2)*storage + (x-alpha(q))/beta(q+1).*b;
  end

  % f = p1*b1 + p0*a0 - p0*b2*beta1/beta2
  f = b/(beta(2)*beta(1)).*(x-alpha(1)) + 1/beta(1)*coeffs(1) - 1/beta(1)*bb*beta(2)/beta(3);
end
