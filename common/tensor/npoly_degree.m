function[P] = npoly_degree(k, d)
% npoly_degree - Returns the number of d-dimensional polynomials of degree k
%
% P = npoly_degree(k,d)
%
%     Returns the number of (independent) d-dimensional polynomials that have
%     total degree k.
%
%     For some reason this is the same as subspace_dimension...wtf.

if k==0
  P = 1;
else
  P = nchoosek(k-1+d, d)*d/k;
end
