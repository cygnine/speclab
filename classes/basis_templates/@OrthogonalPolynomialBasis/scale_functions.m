function p = scale_functions(self,p,n,normalization)
% p = scale_functions(self,p,n,normalization)
%
%     Assumes that the input matrix p has columns that are scaled as
%     L^2-orthornormal evaluations. The vector n indicates which function
%     indices from the basis correspond to the columns of p. This function then
%     rescales each column of p to the FunctionNormalization specification
%     defined by the basis, using the column index ids n.

persistent spdiag orthonormal monic 
if isempty(spdiag)
  from labtools import spdiag
  orthonormal = OrthonormalNormalization.instance();
  monic = MonicNormalization.instance();
end

if not(exist('normalization')==1)
  normalization = self.normalization;
end

if normalization==orthonormal
  % The idea: the affine map induces a Jacobian of map_to_standard_domain.A. If
  % this is already present in the weight rescaling, we don't need to
  % renormalize.
  K = self.scale_weight(1);
  %p = p*sqrt(K/self.map_to_standard_domain.A);
  p = p*sqrt(self.map_to_standard_domain.A/K);
elseif normalization==monic
  N = max(n);
  [a,b] = self.recurrence(0:N);
  b = sqrt(b(:));
  A = [1; repmat(self.map_to_domain.A, [max([N, 0]) 1])];
  b = cumprod(b.*A);
  p = p*spdiag(b(n+1));
else
  error('This basis does not support the given function normalization type');
end
