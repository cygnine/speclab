function[a] = from_naturals(self,n)
% from_naturals -- mapping from naturals for MarginalReverseLexicographicIndexing
% 
% [a] = from_naturals(self,n)
%
%     Computes (0-based) array indexing from (1-based) linear indexing. The
%     output a is a size(n,1) x self.dim matrix, where each row corresponds to the
%     1-based array indexing for the given linear index n.
%
%     Examples:
%          n = [1,2,3,4,5,6,7], self.dim = 2 ----> [0 , 0
%                                                   1 , 1
%                                                   1 , 0
%                                                   0 , 1
%                                                   2 , 2
%                                                   2 , 1
%                                                   2 , 0]
%

persistent subind spdim
if isempty(spdim)
  from speclab.common.tensor import marginalpoly_space_dim as spdim
  from speclab.common.tensor import marginalpoly_subspace_indices as subind
end

n = n(:);
% First sort n:
[n, order] = sort(n);

n_max = n(end);
marginal_degrees = -1*ones(size(n));

% Here's the plan of action:
% 1.) Compute + store all necessary marginal degrees
% 2.) For each marginal degree, compute all multi-indices, pluck out ones in n

% Compute and store marginal degrees:
current_index = 1;
need_higher_degree = true;
N = -1;
while need_higher_degree
  N = N + 1;
  space_dim = spdim(self.dim, N);
  new_index = find(n(current_index:end) > space_dim, 1, 'first');

  % Take care of empty find and reindex to global n-locations:
  if isempty(new_index) 
    % Then all the remaining indices 
    new_index = length(n);
    need_higher_degree = false;
  else
    new_index = current_index + new_index - 2;
  end

  marginal_degrees(current_index:new_index) = N;
  current_index = new_index + 1;
end

% Find markers for new degrees:
degree_inds = [1; find(diff(marginal_degrees)>0)+1; length(n)+1];

% Finally, allocate and create array
a = zeros([numel(n) self.dim]);

for q = 1:length(degree_inds) - 1;
  degree = marginal_degrees(degree_inds(q));
  inds = flipud(subind(self.dim, degree).');

  space_dim = spdim(self.dim, degree-1);

  row_start = degree_inds(q);
  row_end = degree_inds(q+1) - 1;
  a(row_start:row_end,:) = inds(n(row_start:row_end)-space_dim, :);
end

a(order,:) = a;
%a = a.';
a = a + 1;
