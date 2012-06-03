function p = scale_functions(self,p,n,normalization)
% p = scale_functions(self,p,n,normalization)
%
%     Assumes that the input matrix p has columns that are scaled as
%     L^2-orthornormal evaluations. The vector n indicates which function
%     indices from the basis correspond to the columns of p. This function then
%     rescales each column of p to the FunctionNormalization specification
%     defined by the basis, using the column index ids n.

persistent spdiag
if isempty(spdiag)
  from labtools import spdiag
end

if not(exist('normalization')==1);
  normalization = self.normalization;
end

if normalization=='classical'
  % For Discrete Chebyshev polynomials, we'll define this normalization to be the one such
  % that they are normalized on the standard interval with the classical weight
  % to have norm-squared equal to (self.N + n)!/(2*n+1)/(self.N - n - 1)!
  % regardless of affine map.
  temp = gammaln(self.N+n+1) - gammaln(2*n+2) - gammaln(self.N - n);

  temp = gammaln(n+self.alpha+1) + gammaln(n+self.beta+1) - gammaln(n+1) - ...
         gammaln(n+self.alpha+self.beta+1);
  factors = exp(temp);

  p = p*spdiag(sqrt(factors));
else
  p = scale_functions@OrthogonalPolynomialBasis(self, p, n, normalization);
end
