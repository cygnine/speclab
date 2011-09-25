function[x,indices,locations] = nested_filter(x)
% nested_filter -- Removes duplicate nodes for a nested set of nodes
%
% [x, indices, locations] = nested_filter(x)
%   x is an (N x d) matrix representing N points in d-dimensional Euclidean
%   space.  This function removes duplicate nodal entries. If there are K
%   duplicates, the output x is (N-K x d), and indices is a length-N boolean
%   indexing array such that output_x = input_x(indices,:). The third output
%   locations is a list of linear indices of y that the elements in x get mapped
%   to. E.g. to add together weights in a quadrature rule (x,w):
%
%       [y, indices, locations] = nested_filter(x);
%       wy = zeros(size(y));
%       for q = 1:size(y,1)
%         wy(q) = sum(w(locations==q));
%       end


tolerance = 1e-10;

% The metric for coincidence is the l^2 vector norm squared of the difference.
[N,d] = size(x);

if N==1
  indices = true;
  return
end

indices = false([N 1]);
locations = zeros([N 1]);
temp_array = zeros(size(x));
unique_count = 0;
for q = 1:N
  if not(indices(q))
    unique_count = unique_count + 1;

    temp_array(q+1:end,:) = repmat(x(q,:), [N-q 1]) - x(q+1:end,:);
    norms = sqrt(sum(temp_array(q+1:end,:).^2,2));
    indices_matching_current = norms<=tolerance;
    %indices_matching_current = or(indices(q+1:end), norms<=tolerance);
    indices(q+1:end) = or(indices(q+1:end), indices_matching_current);
    locations(q) = unique_count;
    locations([false([q 1]); indices_matching_current]) = q;
  end
end
indices = not(indices);

x = x(indices,:);
