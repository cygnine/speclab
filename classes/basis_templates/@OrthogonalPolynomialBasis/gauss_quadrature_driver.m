function[x,w] = gauss_quadrature_driver(a,b,N)
% gauss_quadrature_driver -- Computes Gauss quadrature from recurrence coefficients
%
% [x,w] = gauss_quadrature_driver(a,b)
%
%     Computes N Gauss quadrature nodes (x) and weights (w) from length-N arrays
%     a and b containing the standard orthonormal recurrence coefficients of the
%     family. This is a static method and can be called without instantiating
%     the class.

a = a(:); b = b(:);
if nargin < 3
  N = min(length(a), length(b));
end

% Form Jacobi matrix
J = spdiags([sqrt([b(2:N);0]) a(1:N) sqrt(b(1:N))], -1:1, N, N);
x = eig(J);
w = OrthogonalPolynomialBasis.evaluate_driver(x, a, b, 0:(N-1), 0);
w = 1./sum(w.^2,2);
