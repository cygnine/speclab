function p = scale_functions(self,p,n,normalization)
% p = scale_functions(self,p,n,normalization)
%
%     Assumes that the input matrix p has columns that are scaled as
%     L^2-orthornormal evaluations. The vector n indicates which function
%     indices from the basis correspond to the columns of p. This function then
%     rescales each column of p to the FunctionNormalization specification
%     defined by the basis, using the column index ids n.

if not(exist('normalization')==1);
  normalization = self.normalization;
end

p = self.scale_functions@MultivariateOrthogonalPolynomialBasis(p, n, normalization);

% Unlike the univariate orthogonal polynomial case, we haven't normalized the
% recurrence coefficient if the user rescales the weight, so we'll have to
% manually rescale things if we normalize functions in a way that depends on the weight.
if normalization=='normal'
  p = p/sqrt(self.scale_weight(1));
end
