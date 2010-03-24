function[coeffs] = exp_gegenbauer_derivatives(m,n,varargin)
% exp_gegenbauer_derivatives -- Evaluates derivatives of expansion coefficients
%
% coeffs = exp_gegenbauer_derivatives(m,n,{lambda=0})
%
%     Assume the following spectral expansion:
%
%          exp(theta*r) = \sum_{n=0}^\infty p_n(r) * \eta_n(theta)
%
%     The polynomials p_n(r) are the L^2-orthonormal Gegenbauer polynomials of
%     order n evaluated at r, and the \eta_n(theta) are the expansion
%     coefficients that depend on the value of theta. (They are evaluated by
%     exp_gegenbauer_expansion.) This function computes D^m \eta_n(0) for inputs
%     m and n and parameter lambda>-1/2. 
%
%     The code is vectorized in m and n, with size(coeffs) = [length(m)
%     length(n)].

persistent strict_inputs spdiag
if isempty(strict_inputs)
  from labtools import strict_inputs spdiag
end

opt = strict_inputs({'lambda'}, {0}, [], varargin{:});
m = m(:);
n = n(:);

assert(all(n>=0), 'Error: the coefficient indices n must be non-negative');
assert(all(m>=0), 'Error: the coefficient indices m must be non-negative');
mmat = repmat(m, [1 length(n)]) ;
nmat = repmat(n.', [length(m) 1]);
mn_valid = (mod(mmat-nmat, 2)==0) & (mmat-nmat)>=0;

n_factors = exp(1/2*log(n+opt.lambda+1/2) + ...
                1/2*gammaln(n+2*opt.lambda+1) - ...
                1/2*gammaln(n+1));

m_factors = sqrt(pi)./(2.^(m+opt.lambda));

coeffs = zeros(size(mmat));
coeffs(mn_valid) = exp(gammaln(mmat(mn_valid)+1) - ...
                       gammaln((mmat(mn_valid) - nmat(mn_valid))/2+1) - ...
                       gammaln((mmat(mn_valid) + nmat(mn_valid))/2 + opt.lambda+3/2));

coeffs = spdiag(m_factors)*coeffs*spdiag(n_factors);
