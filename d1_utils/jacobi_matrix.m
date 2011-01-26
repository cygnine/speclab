function[J] = jacobi_matrix(a, b)
% jacobi_matrix -- Returns the Jacobi matrix
%
% J = jacobi_matrix(a, b)
%
%     Returns the N x N (sparse) Jacobi matrix associated with the polynomial
%     family specified by the length-N recurrence coefficient vectors a and b.

a = a(:); b = b(:);
N = min([length(a), length(b)]);
a = a(1:N);
b = b(1:N);

J = spdiags([sqrt([b(2:end);0]) a sqrt(b)], -1:1, N, N);
