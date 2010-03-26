function[R] = jacobi_rmatrix(N,a,b,A,B);

% [R] = jacobi_rmatrix(N,a,b,A,B);
%
% Creates the full Jacobi R matrix: transporting the modes of (a,b) to the modes
% of (a+A,b+B). A and B must be natural numbers. R is N x N
%
% 20080812 -- acn

A = round(A);
B = round(B);

if (A<0)|(B<0);
  error('Inputs A and B must be non-negative integers');
end

R = speye(N);

for q = 1:A;
  R = jacobi_rsubmatrix(N,a+q,b,'a')*R;
end
for q = 1:B;
  R = jacobi_rsubmatrix(N,a+A,b+q,'b')*R;
end

%tol = 1e-12;
%R(abs(R)<tol) = 0;
R = sparse(R);
