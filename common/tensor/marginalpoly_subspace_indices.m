function[a] = marginalpoly_subspace_indices(dim, k)
% marginalpoly_subspace_indices -- Exponent indices for a polynomial subspace
%
% inds = marginalpoly_subspace_indices(dim, k)
%
%     Returns the exponent indices for the dim-variate polynomial subspace of
%     marginal degree exactly equal to k. Each column of the output is a
%     length-dim exponent multi-index. The indices are ordered
%     lexicographically.
%
%     Example:
%
%                                            
%      marginalpoly_subspace_indices(3,2) ---> 
%                                            
%           

persistent subdim spind subind
if isempty(subdim)
  from speclab.common.tensor import marginalpoly_subspace_dim as subdim
  from speclab.common.tensor import marginalpoly_space_indices as spind
  from speclab.common.tensor import marginalpoly_subspace_indices as subind
end

if (numel(k) > 1) | (numel(dim) > 1)
  error('Support for scalar inputs only');
end


if dim==0
  a = zeros([dim subdim(dim,k)]);
  return
elseif dim==1
  a = k(:);
  return
else
  % Construction actually proceeds by making each multiindex a row vector.
  a = zeros([dim subdim(dim,k)]).';
end

if k > 0 
  row_start = 1;
  row_end = (k+1)^(dim-1);
  a(row_start:row_end,1) = k;
  a(row_start:row_end,2:dim) = flipud(spind(dim-1, k).');

  % Now the next 'set' of rows has k-1 in the first column. Once we generate
  % this, we can just copy+paste for the remainder of the 'sets'. We could do
  % nested function calls, but that makes me...uneasy.
  N = subdim(dim-1,k);
  %submat = zeros([subdim(dim-1,k)] dim-1);
  submat = fliplr(subind(dim-1,k)).';
  for q = (k-1):-1:0
    row_start = row_end+1;
    row_end = row_start + N - 1;
    a(row_start:row_end,1) = q;
    a(row_start:row_end,2:end) = submat;
  end

%  % Column q is of degree k:
%  for q = 2:dim
%    % Number of loops:
%    prefix_inds = flipud(spind(q-1, k-1).');
%    suffix_inds = flipud(spind(dim-q, k).');
%    nloops = k^(q-1); % also is size(prefix_inds,1)
%
%    % For each prefix, copy+paste suffixes
%    for qq = 1:nloops
%      row_start = row_end+1;
%      row_end = row_start + (k+1)^(dim-q) - 1;
%      a(row_start:row_end,1:q-1) = repmat(prefix_inds(qq,:), [(k+1)^(dim-q) 1]);
%      a(row_start:row_end,q) = k;
%      a(row_start:row_end,q+1:dim) = suffix_inds;
%    end
%  end
end

% Things are actually in reverse lexicographic order, plus indices are rows
a = fliplr(a.');
