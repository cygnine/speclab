function[M] = opoly_ip(d, k)
% opoly_ip -- Returns a diagonal matrix
%
% M = opoly_ip(d,k)
%
%     Returns the Euclidean scalar product for the vector space of d-variate
%     polynomials of degree exactly equal k
%
%         ip = v1'*M*v2
%
%     The size of M is determined by the size of the d-variate polynomial space
%     with total degree k.

persistent subdim
if isempty(subdim)
  from speclab.common.tensor import subspace_dimension as subdim
end

M = speye(subdim(d,k));
