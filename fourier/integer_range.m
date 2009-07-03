function[ks] = integer_range(N)
% [KS] = INTEGER_RANGE(N)
%
%     Uses negative-biased integer modal assignment to return a length-N vector
%     containing integer indices. 

if mod(N,2)==0
  ks = -(N/2):((N-2)/2);
else
  ks = -((N-1)/2):((N-1)/2);
end

ks = ks.';
