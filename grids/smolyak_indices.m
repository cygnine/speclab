function[indices] = smolyak_indices(d, k)
% smolyak_indices -- Generates multi-indices required in Smolyak sparse-grid construction
%
% indices = smolyak_indices(d,k)
%
%     Returns an N x d matrix `indices', where each row is a d-dimensional
%     multi-index, where each index is strictly greater than 0. All the
%     multi-indices with total degree k are generated.
%
%     Note therefore that k cannot be less than d (otherwise an empty array is
%     returned).

persistent indexing subdim dim
if isempty(indexing)
  from speclab.common.tensor import linear_to_array_indexing as indexing
  %from speclab.common.tensor import space_dimension as dim
  from speclab.common.tensor import polynomial_space_dimension as dim
  %from speclab.common.tensor import subspace_dimension as subdim
  from speclab.common.tensor import polynomial_subspace_dimension as subdim
end

if k<d
  indices = zeros([0 d]);
  return;
end

N = dim(d,k-d-1);
n = subdim(d,k-d);

indices = indexing(N:(N+n-1), 'dim', d) + 1;
