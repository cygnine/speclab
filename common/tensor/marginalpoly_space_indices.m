function[a] = marginalpoly_space_indices(dim, k)
% marginalpoly_space_indices -- Exponent indices for a polynomial subspace
%
% inds = marginalpoly_space_indices(dim, k)
%
%     Returns the exponent indices for the dim-variate polynomial subspace of
%     marginal degree less than or equal to k. Each column of the output is a
%     length-dim exponent multi-index. The indices are ordered
%     lexicographically.
%
%     Since the size of the resulting matrix is (k+1)^dim x dim, this function
%     quickly gobbles up memory for large dim and k. The workhorse for this
%     function is Matlab's builtin ndgrid function.

if (numel(k) > 1) | (numel(dim) > 1)
  error('Support for scalar inputs only');
end

if dim < 1
  a = zeros([k+1 0]);
  return
end

a = zeros([(k+1)^dim dim]);

outs = cell([dim 1]);
ins = cell([dim 1]);
for q = 1:dim
  ins{q} = 0:k;
end

if dim > 1
  [outs{1:dim}] = ndgrid(ins{:});
else
  outs{q} = (0:k).';
end

for q = 1:dim
  a(:,q) = outs{q}(:);
end
a = fliplr(a).';
