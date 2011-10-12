function[D] = totalpoly_subspace_dim(d,k)
% totalpoly_subspace_dim -- Returns the dimension of a polynomial subspace
%
% D = totalpoly_subspace_dim(d,k)
%
%     Returns the dimension of the space of d-variate polynomials of total
%     degree equal to k.

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
