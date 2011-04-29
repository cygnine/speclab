function[alpha,beta] = modified_chebyshev_algorithm(a, b, m)
% modified_chebyshev_algorithm -- generates recurrence coefficients
%
% [alpha,beta] = modified_chebyshev_algorithm(a, b, m)
%
%     Given the first 2*n recurrence coefficients (a,b) for the polynomial
%     family p, this function generates the first n recurrence coefficients for
%     the polynomial family pi. The length-(2*n) vector m contains integrals of
%     the first (2*n) p polynomials under the measure that defines pi.

N = length(a)/2;

sigmas0 = zeros([2*N 1]);
sigmas1 = reshape(m, [2*N 1]);

alpha = zeros([N 1]);
beta = zeros([N 1]);

alpha(1) = a(1) + m(2)/m(1);
beta(1) = m(1);

toupdate = 2:(2*N-1);

for k = 1:(N-1)
  sigmatemp = sigmas1;
  sigmas1(toupdate) = sigmas1(toupdate+1) - ...
                      (alpha(k) - a(toupdate)).*sigmas1(toupdate) - ...
                      beta(k)*sigmas0(toupdate) + ...
                      b(toupdate).*sigmas1(toupdate-1);
  sigmas0 = sigmatemp;

  alpha(k+1) = a(k+1) + sigmas1(k+2)/sigmas1(k+1) - sigmas0(k+1)/sigmas0(k);
  beta(k+1) = sigmas1(k+1)/sigmas0(k);
  toupdate = (toupdate(1)+1):(toupdate(end)-1);
end
