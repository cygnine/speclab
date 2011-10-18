function[D] = subspace_dimension(d,k)
% subspace_dimension -- Returns the linear dimension of a polynomial space
%
% D = subspace_dimension(d,k)
%
%     Returns the dimension of the space of d-variate polynomials of total
%     degree equal to k.

warning('Use speclab.common.tensor.polynomial_subspace_dimension');

D = zeros(size(k));

for q = 1:numel(k)
  if k(q)==0
    D(q) = 1;
  elseif k(q)<0
    D(q) = 0;
  else
    D(q) = round(d./k(q).*nchoosek(d+k(q)-1, d));
  end
end
