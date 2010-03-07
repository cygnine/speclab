function[D] = space_dimension(d,k)
% space_dimension -- Returns the linear dimension of a polynomial space
%
% D = space_dimension(d,k)
%
%     Returns the dimension of the space of d-variate polynomials of total
%     degree not more than k.

if k<0
  D = 0;
else
  D = nchoosek(d+k,d);
end
