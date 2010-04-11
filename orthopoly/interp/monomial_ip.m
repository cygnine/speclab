function[M] = monomial_ip(d, k)
% monomial_ip -- Returns the de Boor scalar product for monomials
%
% M = monomial_ip(d,k)
%
%     Returns the de Boor scalar product for d-variate multimonomials of degree
%     equal to k. The matrix M that is returned is the symmetric
%     positive-definite matrix M that acts as the bilinear functional taking
%     vectors v1 and v2 to the scalar product:
%
%         ip = v1'*M*v2
%
%     The size of M is determined by the size of the d-variate polynomial space
%     with total degree k.

persistent subdim dim indexing spdiag
if isempty(subdim)
  from labtools import spdiag
  from speclab.common.tensor import space_dimension as dim
  from speclab.common.tensor import subspace_dimension as subdim
  from speclab.common.tensor import linear_to_array_indexing as indexing
end

Msize = subdim(d,k);
%M = zeros(Msize);

% This if clause generates the multidimensional indices for the relevant
% multimonomials. We need them to take factorials.
if k==0
  indices = zeros([1 d]);
else
  d0 = dim(d, k-1);
  indices = indexing(d0:(d0+Msize-1), 'dim', d);
end

%M = spdiag(factorial(k)./prod(factorial(indices), 2));
M = spdiag(1./prod(factorial(indices), 2));
