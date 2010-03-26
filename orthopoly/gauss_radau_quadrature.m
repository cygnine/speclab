function[x,w] = gauss_radau_quadrature(a,b,N,r);
% [x,w] = gauss_radau_quadrature(a,b,N,r);
% 
%     Returns the N-point Gauss-Radau polynomial quadrature for the orthogonal
%     polynomials corresponding to the recurrence coefficients A and B.
%     The fixed Radau point is at the location r.

persistent gq eval_polynomial
if isempty(gq)
  from speclab.orthopoly import gauss_quadrature as gq
  from speclab.orthopoly import eval_polynomial
end

a = a(1:N);
b = b(1:N);

temp = eval_polynomial(r,a,b,[N-2,N-1],'normalization','monic');
a(N) = r - b(N)*temp(1)/temp(2);

a = a(:);
b = b(:);

[x,w] = gq(a,b,N);
