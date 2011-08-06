function[C] = inv_monomial_connection(self, N);
% inv_monomial_connection -- Computes the (inverse) monomial connection matrix
%
% C = inv_monomial_connection(N)
%
%     Computes the N x N lower triangular connection matrix between a particular
%     class of orthogonal polynomials and monomials.  
%
%     The matrix C satisfies:
%
%     p_{n-1} = \sum_{m=1}^{n} C(n,m) x^{m-1}
%
%     where p_m is the m'th degree orthogonal polynomial from this family.
%
%     TODO: this doesn't do anything with scaling/shifting of domains

persistent spdiag 
if isempty(spdiag)
  from labtools import spdiag
end

C = ones(N);
%if opt.dim>1
%  alphas = indexing(0:(N-1), 'dim', opt.dim);
%else
%  alphas = (0:(N-1)).';
%  if isa(recurrence, 'function_handle') % force cell behavior
%    recurrence = {recurrence};
%  end
%end

maxdeg = N-1;
leading_coeffs = zeros([N 1]);

%[a,b] = recurrence{d}(maxdeg+1);
[a,b] = self.recurrence(0:(maxdeg+1));
b = sqrt(b);

leading_coeffs = cumprod(1./b);
%leading_coeffs = leading_coeffs.*local_leading_coeffs(alphas(:,d)+1);

% Build the connection matrix for this dimension and then distribute it
% accordingly
C = diag(leading_coeffs);

a = a(:).';
b = b(:).';
for row = 2:(maxdeg+1);
  % "side" conditions (i.e. those w/o left/right boundary points)
  if row<3
    C(row,1) = 1/b(row)*(-a(row-1)*C(row-1,1));
  else
    C(row,1) = 1/b(row)*(-a(row-1)*C(row-1,1) - b(row-1)*C(row-2,1));

    % All the rest are `vectorizable'
    cols = 2:(row-1);
    C(row,cols) = C(row-1,cols-1) - a(row-1)*C(row-1,cols) - b(row-1)*C(row-2,cols);
    C(row,cols) = C(row,cols)/b(row);
  end
end

% Now distribute
%C = C.*D(alphas(:,d)+1,alphas(:,d)+1);

temp = self.scale_functions(ones([1 size(C,2)]), 0:(size(C,2)-1));
C = spdiag(temp)*C;

%switch self.normalization
%case OrthonomalNormalization.instance()
%  % do nothing
%case MonicNormalization.instance()
%  C = spdiag(1./leading_coeffs)*C;
%end
