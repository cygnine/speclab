function[a] = from_naturals(self,n)
% from_naturals -- mapping from naturals for DegreeIndexing
% 
% [a] = from_naturals(self,n)
%
%     Computes degree in dim-dimensional multi-index space that corresponds to
%     input linear indexing. 
%
%     Examples:
%          n = [1,2,3], self.dim = 2 ----> [0, 1, 1]

if self.dim==1
  a = n-1;
  return
end

assert(all(n > 0), 'Error: inputs must be strictly positive');

% First find maximum degree
max_degree = 0;
while any(n > self.total_polynomial_space_dimension(self.dim, max_degree));
  max_degree = max_degree + 1;
end

% Now, generate histogram bins:
edges = 1 + self.total_polynomial_space_dimension(self.dim, 0:(max_degree+1));
[garbage, a] = histc(n, edges);
