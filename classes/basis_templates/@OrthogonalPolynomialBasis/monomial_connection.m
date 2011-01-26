function[C] = monomial_connection(self, N);
% monomial_connection -- Computes connection coefficients to monomials
%
% C = monomial_connection(N)
%
%     Computes the N x N lower triangular connection matrix between monomials
%     and a the orthogonal polynomial class.
%
%     The matrix C satisfies:
%
%     x^{n-1} = \sum_{m=1}^{n} C(n,m) p_{m-1},
%
%     where p_m is the m'th degree orthogonal polynomial.
%
%     TODO: this doesn't do anything wtih scaling/shifting of domains yet.

persistent spdiag 
if isempty(spdiag)
  from labtools import spdiag
end

%opt = strict_inputs({'dim', 'normalization', 'weight_normalization'}, {1, 'normal', ''}, [], varargin{:});

C = ones(N);
%if opt.dim>1
%  alphas = indexing(0:(N-1), 'dim', opt.dim);
%else
%  alphas = (0:(N-1)).';
%  if isa(recurrence, 'function_handle') % force cell behavior
%    recurrence = {recurrence};
%  end
%end

ns = (0:(N-1)).';
maxdeg = N-1;

[a,b] = self.recurrence(maxdeg+1);
b = sqrt(b);

leading_coeffs = cumprod(1./b);
%leading_coeffs = leading_coeffs.*local_leading_coeffs(alphas(:,d)+1);

% Build the connection matrix for this dimension and then distribute it
% accordingly
C = zeros(N);
C = diag(1./leading_coeffs);

a = a(:).';
b = b(:).';
for row = 2:(maxdeg+1);
  % "side" conditions (i.e. those w/o left/right boundary points)
  C(row,1) = a(1)*C(row-1,1) + b(2)*C(row-1,2);

  % All the rest are `vectorizable'
  cols = 2:(row-1);
  C(row,cols) = a(cols).*C(row-1,cols) + b(cols).*C(row-1,cols-1) + ...
                b(cols+1).*C(row-1,cols+1);
end

% Now distribute
%C = C.*D(alphas(:,d)+1,alphas(:,d)+1);

switch self.normalization
case OrthonormalNormalization.instance()
  % do nothing
case MonicNormalization.instance()
  C = C*spdiag(1./leading_coeffs);
otherwise
  error('Normalization not yet supported');
end
