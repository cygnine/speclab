function[y] = lebfun(x,z)
% lebfun -- Evaluates the Lebesgue function of a one-dimensional set of nodes
%
% y = lebfun(x,z)
%
%     For a vector x, evaluates the Lebesgue function associated with the nodes
%     x at the locations z.
%
%     This function is speed-optimized, meaning that if length(x) times length(z)
%     is larger than your memory allows, Matlab will throw an error.

x = sort(x(:));
zsize = size(z);
z = z(:);

y = zeros(size(z));

N = length(x);
Nz = length(z);

% (z-x_j)
linear_factors = repmat(z, [1 N]) - repmat(x.', [Nz 1]);

indexing = repmat((1:N).', [N 1]);
diagonal_indices = (0:(N-1))*N + (1:N);
indexing(diagonal_indices) = [];
indexing = reshape(indexing, [N-1 N]);

for j = 1:N;
  const = prod(x(j) - x(indexing(:,j)));
  y = y + abs(prod(linear_factors(:, indexing(:,j)),2)/const);
end

y = reshape(y, zsize);
