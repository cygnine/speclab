function[x,w] = gauss_lobatto_quadrature(a,b,N,r1,r2);
% [X,W] = GAUSS_LOBATTO_QUADRATURE(A,B,N,R1,R2);
%
%     Returns the N-point Gauss-Lobatto polynomial quadrature for the orthogonal
%     polynomials corresponding to the recurrence coefficients a and b.
%     The fixed Lobatto points are r1, r2.

global handles;
gq = handles.speclab.OrthogonalPolynomial1D.gauss_quadrature;

a = a(1:N);
b = b(1:N);

temp = eval_opoly([r1; r2],a,b,[N-1,N-2]);
modif = inv(temp)*[r1*temp(1,1); r2*temp(2,1)];

% Lobatto modification for Jacobi matrix: 1999_gautschi
a(N) = modif(1);
b(N) = modif(2);

a = a(:);
b = b(:);

[x,w] = gq(a,b,N);

%J = spdiags([sqrt([b(2:end);0]) a sqrt(b)], -1:1, N, N);
%
%[v,d] = eig(full(J));
%x = diag(d);
%w = b(1)*v(1,:).^2.';
