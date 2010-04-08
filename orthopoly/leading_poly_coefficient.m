function[k] = leading_poly_coefficient(n, recurrence, dim, normalization, weight_normalization, poly_parameters)
% k = leading_poly_coefficient(n, recurrence, dim, normalization, weight_normalization, poly_parameters)
%
%     Evaluates the leading coefficient of orthogonal polynomials.
%
%     The input 'recurrence' is a function handle pointing to the recurrence
%     relation for the orthogonal polynomial family.
%
%     'poly_parameters' is a struct with fields given by any necessary
%     parameters for the orthogonal polynomial family.
%
%     'normalization' and 'weight_normalization' are strings that determine how
%     the polynomials are normalized ('normalization') with respect to the
%     particular form of the weight_function ('weight_normalization').

persistent indexing
if isempty(indexing)
  from speclab.common.tensor import linear_to_array_indexing as indexing
end

if dim==1
  n_array = n(:); % no need to call linear_to_array_indexing
  N = max(n)+2;
else
  n_array = indexing(n,'dim',dim);
  N = max(max(n_array)) + 1;
end

k = ones([size(n_array, 1) 1]);

switch normalization
case {'normal', 'orthonormal'}
  for q = 1:dim
    [a,b] = recurrence(N, poly_parameters(q));
    b = cumprod(sqrt(b));

    k = k.*b(n_array(:,q)+1);
  end
  k = 1./k;

case 'monic' % super easy
  k = ones(size(n));
otherwise
  error('Normalization type not supported')
end

switch weight_normalization
case ''
  % do nothing
case 'probability'
  error('Not yet implemented');
otherwise
  error('Weight normalization type not supported');
end

k = reshape(k, size(n));
