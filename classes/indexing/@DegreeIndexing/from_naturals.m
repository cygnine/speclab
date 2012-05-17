function[a] = from_naturals(self,n)
% from_naturals -- mapping from naturals for DegreeIndexing
% 
% [a] = from_naturals(self,n)
%
%     Computes (0-based) array indexing from (1-based) linear indexing. The
%     output a is a size(n,1) x self.dim matrix, where each row corresponds to the
%     1-based array indexing for the given linear index n.
%
%     Examples:
%          n = [1,2,3], self.dim = 2 ----> [0 , 0
%                                           1 , 0
%                                           0 , 1
%                                           2 , 0
%                                           1 , 1
%                                           0 , 2]
%

persistent subind spdim subdim
if isempty(spdim)
  from speclab.common.tensor import totalpoly_subspace_dim as subdim
  from speclab.common.tensor import totalpoly_space_dim as spdim
  from speclab.common.tensor import totalpoly_subspace_indices as subind
end

if self.dim==1
  a = n-1;
  return
end

assert(all(n > 0), 'Error: inputs must be strictly positive');

n = n(:);
% First sort n:
[n_sorted, order] = sort(n);
degree_inds = [1; find(diff(n_sorted)>0)+1; length(n_sorted)+1];

% Allocate total array
a_row_sizes = subdim(self.dim, n-1);
a_row_indices = [0; cumsum(a_row_sizes)];
a = zeros([sum(a_row_sizes) self.dim]);
for q = 1:length(degree_inds)-1
  indices = subind(self.dim, n_sorted(degree_inds(q))-1).';

  % For each reptition of the degree, paste into array a
  for qq = 1:(degree_inds(q+1) - degree_inds(q))
    % Find starting index in a_row_indices
    index = order(qq - 1 + degree_inds(q));
    row_inds = (a_row_indices(index)+1):(a_row_indices(index+1));
    a(row_inds,:) = indices;
  end
end
