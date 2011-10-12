function[a] = totalpoly_subspace_indices(dim, k)
% totalpoly_subspace_indices -- Exponent indices for a polynomial subspace
% 
% inds = totalpoly_subspace_indices(dim, k)
%
%     Returns the exponent indices for the dim-variate polynomial subspace of
%     total degree exactly equal to k. Each column of the output is a length-dim
%     exponent multi-index. The indices are ordered lexicographically.
%
%     Example:
%
%                                            [ 0     0     0     1     1     2 ]
%       totalpoly_subspace_indices(3,2) ---> [ 0     1     2     0     1     0 ]
%                                            [ 2     1     0     1     0     0 ]
%           
%            

persistent subdim
if isempty(subdim)
  from speclab.common.tensor import totalpoly_subspace_dim as subdim
end

if (numel(k) > 1) | (numel(dim) > 1)
  error('Support for scalar inputs only');
end

% Construction actually proceeds by making each multiindex a row vector.
a = zeros([dim subdim(dim,k)]).';

if k > 0
  row_id = 1;
  current_row = zeros([1 dim]);
  current_row(1) = k;
  a(row_id,:) = current_row;
  row_id = row_id + 1;

  % "The traveling ones-man method"
  onesman_home = 1;
  onesman_location = 1;

  finished = false;

  while not(finished);
    onesman_pilgrimage();
  end

  % Rows are now ordered reverse-lexicographically.
  a = flipud(a);
end

a = a.';

function[] = onesman_pilgrimage();
  while onesman_location < dim
    onesman_location = onesman_location + 1;
    current_row(onesman_location-1) = current_row(onesman_location-1) - 1;
    current_row(onesman_location) = current_row(onesman_location) + 1;
    a(row_id,:) = current_row;
    row_id = row_id + 1;
  end

  if onesman_home + 1 == dim 
    % Then make all the other onesman in column dim-1 travel as well
    while current_row(onesman_home)>0
      current_row(end) = current_row(end) + 1;
      current_row(end-1) = current_row(end-1) - 1;
      a(row_id,:) = current_row;
      row_id = row_id + 1;
    end
  end

  if current_row(end)==k
    finished = true;
    return % done!
  end

  % Now update new home for (next) onesman
  % There must exist a zero in some column; find the last consecutive one from
  % the right

  columns = find(current_row, 2, 'last');
  current_row(columns(1)) = current_row(columns(1)) - 1;
  current_row(columns(1)+1) = current_row(end) + 1;
  current_row(end) = 0;
  a(row_id,:) = current_row;
  row_id = row_id + 1;

  onesman_home = columns(1)+1;
  onesman_location = columns(1)+1;
end

end
