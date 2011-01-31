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
  % For Chebyshev polynomials, this is the normalization where the L^2 norm
  % squared of P_0 is pi and for P_n is pi/2 for all n>0.
  factors = pi/2*ones(size(n));
  p = p*spdiag(sqrt(factors));
else
  p = self.scale_functions@OrthogonalPolynomialBasis(p, n, normalization);
end
