function[lmbda] = derivative_expansion(self, n, d)
% derivative_expansion -- expansion coefficients of basis derivatives
%
% lmbda = derivative_expansion(self, n, d)
%
%     n is a vector of valid indices, and d is an array, where each row is a
%     multi-index of positive integers.
%
%     If size(d,1)==1: returns a length(n) x M matrix with expansion
%     coefficients of self such that the d'th derivative of the n'th orthogonal
%     polynomials has expansion coefficients lmbda. Row q+1 of lmbda has
%     coefficients:
%
%         p^{(d)}_{n(q)} = \sum_{m=1}^{M} lmbda(q,m+1) * p_{m},
%
%     where for clarity above we have assumed that n and m are
%     OneBasedIndexing. p^{(d)}_n is the d'th (partial) derivative of the basis
%     function p_n associated with self. The number of columns M of lmbda is
%     chosen to be as small as possible: M = max(n+1) - |d|.
%
%     If d is a matrix: returns a three-dimensional array lmbda, where
%     lmbda(:,:,r) is the matrix lmbda corresponding to d = d(r,:). The number of
%     columns of lmbda is the numbers of columns for the smallest value of d.

% First change indices to ZeroBasedIndexing on each dimension:
[n, nsize, numeln] = self.indexing(n);

% Here's the plan:
% - see which dimensions the multi-index d is taking derivatives in
% - compute univariate differentiation matrices for those dimensions
% - pass derivative coefficients into array of tensorized indexing

% First: loop over each dimension, storing the appropriate derivative matrix
% Find max_degree in each dimenion:
max_degree = max(n, [], 1);
% Compute which derivatives are necessary in each dimension:
needed_derivatives = zeros(size(d));
temp = sort(d, 1);
for q = 1:self.dim
  dimension_derivatives = setdiff(unique(temp(:,q)), 0);
  needed_derivatives(1:length(dimension_derivatives),q) = dimension_derivatives;
end
max_derivative = max(needed_derivatives(:));

% Now compute univariate matrices for those derivatives:
univariate_matrices = cell([max_derivative+1 self.dim]);
for q = 1:self.dim
  counter = 1;
  %number_derivatives = find(needed_derivatives(:,q) == 0, 1, 'first') - 1;
  number_derivatives = sum(needed_derivatives(:,q) > 0);

  % Get all the matrices for this dimension in one fell swoop
  temporary_matrix_storage = self.bases{q}.derivative_expansion(0:max_degree(q), needed_derivatives(1:number_derivatives,q));
  % And do awkward looping to distribute them into the cell array
  for qq = 1:number_derivatives
    univariate_matrices{needed_derivatives(qq,q), q} = temporary_matrix_storage(:,:,qq);
  end
  univariate_matrices{max_derivative+1,q} = speye(max_degree(q)+1);
end

% Allocate sparse array for expansion coefficients
lmbda = cell([size(d,1) 1]);
for q = 1:size(d,1);
  % Determine size of matrix for this multi-index
  N = self.internal_indexing.inv(max_degree - d(q,:));
  lmbda{q} = spalloc(size(n,1), N, self.dim*size(n,1));
end

% Now for any element of d that equals zero, we just need the univariate
% connection matrix to be the identity. A not-so-elegant solution to this is to
% set all d's that are 0 to be a higher index than derivative we'll need
% (max_derivative+1). Then we can index univariate_matrices like usual.
d(d==0) = max_derivative+1;

% For each row of n, we'll make a CombinationGenerator to see where to place the
% coefficients. This is *super* slow, but I don't know an easier way to code it
% up.
indices = cell([self.dim 1]); % temporary storage array 
for row = 1:size(n,1);
  % Do stuff
  for multiindex = 1:size(d,1);
    % Find indices necessary in each dimension
    % Also is better to store necessary univariate coeffs in an array to
    % facilitate easy indexing later on
%    coeffs = zeros([
    compute = true;
    for dim = 1:self.dim
      indices{dim} = find(abs(univariate_matrices{d(multiindex,dim), dim}(n(row,dim)+1,:)) > 0) - 1;
      compute = not(isempty(indices{dim}));
    end

    % If the derivative doesn't make the function vanish:
    if compute
      G = CombinationGenerator(indices{:});
      % Each instance from G is a multi-index corresponding to expansion in the basis
      inds = G.next(G.set_cardinality);
      cols = self.internal_indexing.inv(inds); % indexing to 1, 2, 3...

      % Initialize:
      lmbda{multiindex}(row,cols) = 1;
      % Now for each columns of inds, we'll pluck that coefficient from the
      % univariate matrix, and multiplex it into row "row" of lmbda{multiindex}.
      for dim = 1:self.dim;
        lmbda{multiindex}(row,cols) = lmbda{multiindex}(row,cols).*univariate_matrices{d(multiindex,dim),dim}(n(row,dim)+1,inds(:,dim)+1);
      end
    end
  end
end

% Remove some sizes of the matrix
for multiindex = 1:size(d,1);
  temp = find(sum(abs(lmbda{multiindex}),1)>0, 1, 'last');
  lmbda{multiindex}(:,temp+1:end) = [];
end
if size(d,1)==1
  lmbda = lmbda{1};
end
