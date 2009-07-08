function[x,w] = gauss_radau_quadrature(a,b,N,r);
% [X,W] = GAUSS_RADAU_QUADRATURE(A,B,N,R);
% 
%     Returns the N-point Gauss-Radau polynomial quadrature for the orthogonal
%     polynomials corresponding to the recurrence coefficients A and B.
%     The fixed Radau point is at the location r.

global handles;
opoly = handles.speclab.OrthogonalPolynomial1D;
gq = handles.speclab.OrthogonalPolynomial1D.gauss_quadrature;

a = a(1:N);
b = b(1:N);

temp = opoly.eval_polynomial(r,a,b,[N-2,N-1],'normalization','monic');
a(N) = r - b(N)*temp(1)/temp(2);

a = a(:);
b = b(:);

[x,w] = gq(a,b,N);

%J = spdiags([sqrt([b(2:end);0]) a sqrt(b)], -1:1, N, N);
%
%[v,d] = eig(full(J));
%x = diag(d);
%w = b(1)*v(1,:).^2.';
