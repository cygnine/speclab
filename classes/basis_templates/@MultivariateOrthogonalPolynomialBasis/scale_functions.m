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

if not(exist('normalization')==1)
  normalization = self.normalization;
end

if normalization=='orthonormal'
  % We only need to care about affine mapping: by default self.recurrence
  % already has the weight normalization information built in.
  p = p*sqrt(det(self.map_to_standard_domain.A));

  % Below: the old way before information was built into the recurrence relation.
  % The idea: the affine map induces a Jacobian of map_to_standard_domain.A. If
  % this is already present in the weight rescaling, we don't need to
  % renormalize.
  %K = self.scale_weight(1);
  %p = p*sqrt(self.map_to_standard_domain.A/K);
else
  error('This basis does not support the given function normalization type');
end
