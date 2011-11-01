function[D] = polynomial_space_dimension(d,k)
% polynomial_space_dimension -- Returns the linear dimension of a polynomial space
%
% D = polynomial_space_dimension(d,k)
%
%     Returns the dimension of the space of d-variate polynomials of total
%     degree less than or equal to k.

D = zeros(size(k));

for q = 1:numel(k);
  if k(q)<0
    D(q) = 0;
  else
    D(q) = nchoosek(d+k(q),d);
  end
end
