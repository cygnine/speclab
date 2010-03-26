function[f] =jacobi_rmatrix_invert(F,R);

% [f] =jacobi_rmatrix_invert(F,R);
% f = inv(R)*F
% 
% Inverts the R matrix using a back-substitution sequential method. Takes
% advantage of the special structure of the R matrix. (Banded, upper-triangular)
% Is vectorized in F i.e. F can have multiple columns.
%
% 20080813 -- acn

N = size(R,1);

offset = max(find(R(1,:)));
f = zeros(size(F));

f(N,:) = 1/R(N,N)*F(N,:);

for q = (N-1):(-1):(N-offset+1);
  f(q,:) = 1/R(q,q)*(F(q,:) - R(q,(q+1):end)*f((q+1):end,:));
end

for q = (N-offset):(-1):1;
  f(q,:) = 1/R(q,q)*(F(q,:) - R(q,(q+1):(q+offset))*f((q+1):(q+offset),:));
end
