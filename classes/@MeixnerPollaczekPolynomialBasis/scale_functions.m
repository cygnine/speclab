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
  % For Meixner-Pollaczek polynomials, we use the following analytic expression
  % for the L^2 norms.
  factors = 1/(2*sin(self.phi))^(2*self.lambda);
  factors = factors.*exp(gammaln(n + 2*self.lambda) - gammaln(n+1));
  p = p*spdiag(sqrt(factors));
else
  p = self.scale_functions@OrthogonalPolynomialBasis(p, n, normalization);
end
