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

if isa(normalization, 'ProbabilistFunctionNormalization');
  % For Hermite polynomials, these are monic on the standard interval. For
  % mapped intervals we'll defined them as the maps of the monic functions (no
  % rescaling).

  % Copy-paste monic normalization
  p = self.scale_functions@OrthogonalPolynomialBasis(p, n, MonicNormalization.instance());
elseif isa(normalization, 'PhysicistFunctionNormalization') | ...
       isa(normalization, 'ClassicalFunctionNormalization')
  % On the standard interval, these are functions whose leading coefficient is
  % 2^n. I.e. first they're monic, then we scale them by 2^n.
  % On mapped intervals, we just map the functions -- no scaling is done.
  N = max(n);
  [a,b] = self.recurrence(0:N);
  b = sqrt(b(:));
  A = [1; repmat(self.map_to_domain.A, [max([N, 0]) 1])];
  b = cumprod(b.*A);
  p = p*spdiag(b(n+1));
  p = p*spdiag(2.^(n));
end
