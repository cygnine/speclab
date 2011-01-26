function[J] = jacobi_matrix(self, N)
% jacobi_matrix -- Returns the Jacobi matrix
%
% J = jacobi_matrix(self, N)
%
%     Returns the N x N (sparse) Jacobi matrix associated with the polynomial
%     family.

persistent jmat
if isempty(jmat)
  from speclab.d1_utils import jacobi_matrix as jmat
end

[a,b] = self.recurrence(0:(N-1));
J = jmat(a, b);

%a = a(:); b = b(:);
%J = spdiags([sqrt([b(2:end);0]) a sqrt(b)], -1:1, N, N);
