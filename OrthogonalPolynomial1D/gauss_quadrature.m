function[x,w] = gauss_quadrature(a,b,N);
% [X,W] = GAUSS_QUADRATURE(A,B,N);
%
%     Returns the N-point Gaussian polynomial quadrature for the orthogonal
%     polynomials corresponding to the recurrence coefficients A and B.

a = a(1:N);
b = b(1:N);

a = a(:);
b = b(:);

J = spdiags([sqrt([b(2:end);0]) a sqrt(b)], -1:1, N, N);

[v,d] = eig(full(J));
x = diag(d);
w = b(1)*v(1,:).^2.';
