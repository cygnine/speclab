function[x,w] = gauss_quadrature(a,b);
% [x,w] = gauss_quadrature(a,b);
%
%     Returns the N-point Gauss polynomial quadrature for the orthogonal
%     polynomials corresponding to the recurrence coefficients given in the
%     length-N vectors a and b.

persistent jacobi_matrix opoly_evaluate
if isempty(jacobi_matrix)
  from speclab.d1_utils import jacobi_matrix opoly_evaluate
end

J = jacobi_matrix(a, b);
n = size(J,1);
x = eig(J);
w = 1./sum(opoly_evaluate(x, a, b, 0:(n-1), 0).^2, 2);
