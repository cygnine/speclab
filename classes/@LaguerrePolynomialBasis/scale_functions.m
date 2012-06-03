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
  % For Laguerre polynomials, these are normalized so that their squared L^2
  % norm is equal to gamma(n+alpha+1)/gamma(n+1).

  factors = exp(gammaln(n + self.alpha + 1) - gammaln(n + 1));
  p = full(p*spdiag(sqrt(factors)));
else
  p = scale_functions@OrthogonalPolynomialBasis(self, p, n, normalization);
end
