function[ks] = integer_range(N)
% [ks] = integer_range(N)
%
%     Creates the range of integers from -N/2, -N/2+1, ..., N/2, with a negative
%     bias if N is even. 

if mod(N,2)==0
  negative_limit = -N/2;
  positive_limit = (N-2)/2;
else
  negative_limit = -(N-1)/2;
  positive_limit = (N-1)/2;
end

ks = [negative_limit:positive_limit].';
