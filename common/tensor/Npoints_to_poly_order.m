function[k] = Npoints_to_poly_order(d, N)
% Npoints_to_poly_order -- Returns the minimum polynomial degree for N points
%
% k = Npoints_to_poly_order(d,N)
%
%     In d-dimensional space, this function returns the total degree k
%     polynomial space necessary so that dim(degree k or less polys) >= N.

persistent dim
if isempty(dim)
  from speclab.common.tensor import space_dimension as dim
end

k = 0;
while dim(d, k)<N
  k = k + 1;
end
