function[C] = symmetric_jacobi_ip(d,k, varargin)
% symmetric_jacobi_ip -- The symmetric Jacobi polynomial inner product
%
% C = symmetric_jacobi_ip(d,k, {lambda=0})
%
%     Returns a D(d,k) x 1 vector that define the (scaled) least polynomial
%     inner product using orthonormal symmetric Jacobi polynomials of class
%     (lambda,lambda). D(d,k) is the size of the (exact) degree-k polynomial
%     space in d dimensions. 
%
%     To evaluate the inner product between two vectors v and w then is simply
%     sum(v.*w.*C).

persistent coeffs strict_inputs
persistent sp_dim subsp_dim indexing monomial_ip
if isempty(strict_inputs)
  from labtools import strict_inputs
  from speclab.apps import exp_gegenbauer_derivatives as coeffs
  from speclab.common.tensor import subspace_dimension as subsp_dim
  from speclab.common.tensor import space_dimension as sp_dim
  from speclab.common.tensor import linear_to_array_indexing as indexing
  from speclab.orthopoly.interp import monomial_ip
end

opt = strict_inputs({'lambda'}, {zeros([d 1])}, [], varargin{:});

N = subsp_dim(d,k);
C = ones([N 1]);

alphas = indexing(sp_dim(d,k-1):sp_dim(d,k)-1, 'dim', d);

for q= 1:d
  cfs = diag(coeffs(0:k,0:k, 'lambda', opt.lambda(q)));

  C = C.*cfs(alphas(:,q)+1);
end

C = C(1,1)./C;  % Normalize the ip so C(1) = 1
M = diag(monomial_ip(d,k));  % This is already normalized so M(1) = 1

C = C.^2.*M;
