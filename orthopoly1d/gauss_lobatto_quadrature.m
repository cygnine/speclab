function[x,w] = gauss_lobatto_quadrature(a,b,N,r1,r2);
% [x,w] = gauss_lobatto_quadrature(a,b,N,r1,r2);
%
%     Returns the N-point Gauss-Lobatto polynomial quadrature for the orthogonal
%     polynomials corresponding to the recurrence coefficients a and b.
%     The fixed Lobatto points are r1, r2.

persistent gq eval_polynomial
if isempty(gq)
  from speclab.orthopoly1d import gauss_quadrature as gq
  from speclab.orthopoly1d import eval_polynomial
end

a = a(1:N);
b = b(1:N);

temp = eval_polynomial([r1; r2],a,b,[N-1,N-2],'normalization','monic');
modif = inv(temp)*[r1*temp(1,1); r2*temp(2,1)];

% Lobatto modification for Jacobi matrix: 1999_gautschi
a(N) = modif(1);
b(N) = modif(2);

a = a(:);
b = b(:);

[x,w] = gq(a,b,N);
