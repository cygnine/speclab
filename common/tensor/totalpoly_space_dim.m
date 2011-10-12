function[D] = totalpoly_space_dim(d,k)
% totalpoly_space_dim -- Returns the linear dimension of a polynomial space
%
% D = totalpoly_space_dim(d,k)
%
%     Returns the dimension of the space of d-variate polynomials of total
%     degree not more than k.

D = zeros(size(k));

for q = 1:numel(k);
  if k(q)<0
    D(q) = 0;
  else
    D(q) = nchoosek(d+k(q),d);
  end
end
