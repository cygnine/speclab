function[x,indices] = nested_filter(x)
% nested_filter -- Removes duplicate nodes for a nested set of nodes
%
% [x, indices] = nested_filter(x)
%   x is an (N x d) matrix representing N points in d-dimensional Euclidean
%   space.  This function removes duplicate nodal entries. If there are K
%   duplicates, the output x is (N-K x d), and indices is a length-N boolean
%   indexing array such that output_x = input_x(indices,:).

tolerance = 1e-10;

% The metric for coincidence is the l^2 vector norm squared of the difference.
[N,d] = size(x);

if N==1
  indices = true;
  return
end

N_removed_nodes = 0;
indices = false([N 1]);
temp_array = zeros(size(x));
for q = 1:N
  if not(indices(q))
    temp_array(q+1:end,:) = repmat(x(q,:), [N-q 1]) - x(q+1:end,:);
    norms = sqrt(sum(temp_array(q+1:end,:).^2,2));
    indices(q+1:end) = or(indices(q+1:end), norms<=tolerance);
  end
end
indices = not(indices);

x = x(indices,:);
