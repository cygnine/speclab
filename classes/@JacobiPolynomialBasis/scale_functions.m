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

if isa(normalization, 'ClassicalFunctionNormalization');
  % For Jacobi polynomials, we'll define this normalization to be the one such
  % that p_n(R) = nchoosek(n+alpha, n), where R is the right-hand endpoint,
  % regardless of affine map.
  temp = gammaln(n+self.alpha+1) + gammaln(n+self.beta+1) - gammaln(n+1) - ...
         gammaln(n+self.alpha+self.beta+1);
  factors = 2^(self.alpha+self.beta+1)./(2*n+self.alpha+self.beta+1).*exp(temp);
  p = p*spdiag(sqrt(factors));
else
  p = self.scale_functions@OrthogonalPolynomialBasis(p, n, normalization);
end
