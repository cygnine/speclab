function[C] = jacobi_monomial_connection(k, varargin)
% jacobi_monomial_connection -- The multidimensional Jacobi-monomial connection matrix
%
% C = jacobi_monomial_connection(k,{dim=1, lambda=zeros(size(dim))})
%
%     Computes the D(dim,k) x D(dim,k) connection matrix between tensor-product
%     Jacobi polynomials and monomials. D(dim,k) is the size of the space
%     spanned by k-th degree polynomials (or less) in dim dimensions.  The
%     Jacobi polynomials in dimension q are symmetric of class (lambda(q),
%     lambda(q)).
%
%     Given a vector v of Jacobi polynomial coefficients, the corresponding
%     monomial coefficients are given by C*v.

persistent coeffs strict_inputs
persistent sp_dim subsp_dim indexing
if isempty(strict_inputs)
  from labtools import strict_inputs
  from speclab.apps import exp_gegenbauer_derivatives as coeffs
  from speclab.common.tensor import space_dimension as sp_dim
  from speclab.common.tensor import linear_to_array_indexing as indexing
end

opt = strict_inputs({'dim', 'lambda'}, {1, 0}, [], varargin{:});
if length(opt.lambda)>opt.dim
  opt.lambda = ones([opt.dim 1])*opt.lambda;
end

N = sp_dim(opt.dim,k);
C = ones(N);
alphas = indexing(0:(N-1), 'dim', opt.dim);

for q = 1:opt.dim
  C_dim1 = coeffs(0:k, 0:k, 'lambda', opt.lambda(q)).';

  C = C.*C_dim1(alphas(:,q)+1,alphas(:,q)+1);
end
