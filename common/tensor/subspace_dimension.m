function[D] = subspace_dimension(d,k)
% subspace_dimension -- Returns the linear dimension of a polynomial space
%
% D = subspace_dimension(d,k)
%
%     Returns the dimension of the space of d-variate polynomials of total
%     degree equal to k.

if k==0
  D = 1;
elseif k<0
  D = 0;
else
  D = round(d./k.*nchoosek(d+k-1, d));
end
