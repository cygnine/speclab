function[C] = monomial_connection_driver(N, recurrence, varargin)
% monomial_connection_driver -- Computes the monomial connection matrix
%
% C = monomial_connection_driver(N, recurrence, {dim=1, normalization='normal', weight_normalization=''})
%
%     Computes the N x N lower triangular connection matrix between monomials
%     and a particular class of orthogonal polynomials. The input `recurrence'
%     is a function handle with the syntax [a,b] = recurrence(N) generating the
%     first N recurrence coefficients for an orthogonal polynomial sequence.
%
%     The matrix C satisfies:
%
%     x^{n-1} = \sum_{m=1}^{n} C(n,m) p_{m-1},
%
%     where p_m is the m'th degree orthogonal polynomial that is generated
%     by the recurrence relation. The optional inputs `normalization' and
%     `weight_normalization' specify how the polynomials are normalized.
%
%     The optional input `dim' specifies the dimension. The input `recurrence'
%     should be a length-dim cell array of function handles if dim>1.

persistent strict_inputs spdiag indexing
if isempty(strict_inputs)
  from labtools import strict_inputs spdiag
  from speclab.common.tensor import linear_to_array_indexing as indexing
end

opt = strict_inputs({'dim', 'normalization', 'weight_normalization'}, {1, 'normal', ''}, [], varargin{:});

C = ones(N);
if opt.dim>1
  alphas = indexing(0:(N-1), 'dim', opt.dim);
else
  alphas = (0:(N-1)).';
  if isa(recurrence, 'function_handle') % force cell behavior
    recurrence = {recurrence};
  end
end

maxdeg = max(max(alphas));
leading_coeffs = zeros([N 1]);

for d = 1:opt.dim
  [a,b] = recurrence{d}(maxdeg+1);
  b = sqrt(b);

  local_leading_coeffs = cumprod(1./b);
  leading_coeffs = leading_coeffs.*local_leading_coeffs(alphas(:,d)+1);
  % Build the connection matrix for this dimension and then distribute it
  % accordingly
  D = zeros(maxdeg+1);
  D = diag(1./local_leading_coeffs);

  a = a.';
  b = b.';
  for row = 2:(maxdeg+1);
    % "side" conditions (i.e. those w/o left/right boundary points)
    D(row,1) = a(1)*D(row-1,1) + b(2)*D(row-1,2);

    % All the rest are `vectorizable'
    cols = 2:(row-1);
    D(row,cols) = a(cols).*D(row-1,cols) + b(cols).*D(row-1,cols-1) + ...
                  b(cols+1).*D(row-1,cols+1);
  end

  % Now distribute
  C = C.*D(alphas(:,d)+1,alphas(:,d)+1);
end

switch opt.normalization
case {'normal', 'orthonormal'}
  % do nothing
case 'monic'
  C = C*spdiag(1./leading_coeffs);
end
