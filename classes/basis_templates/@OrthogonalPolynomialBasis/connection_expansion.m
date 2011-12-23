function[v] = connection_expansion(self, u, other)
% connection_expansion -- Rewrites a polynomial expansion
%
% v = connection_expansion(self, u, other)
%
%     u is an array with N rows, where each column contains expansion
%     coefficients in the basis "other". This function translates each column
%     into expansion coefficients in the basis "self". 

if nargin < 3
  v = u;
  return;
end

N = size(u,1);
% All we need are the connection coefficients: 
C = other.orthogonal_connection(self,N);

v = C'*u;
