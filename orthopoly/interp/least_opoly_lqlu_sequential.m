function[l,u,p,v,k] = least_opoly_lqlu_sequential(theta, varargin)
% tensor_lu -- Computes de Boor's LU factorization with LQ factorizations
%
% [l,u,v,p,k] = least_opoly_lqlu_sequential(theta, [tol=1e-10, ip=eye(),  basis='legendre'])
%
%     Does the same thing as least_opoly_lqlu except performs computations one
%     row at a time.

persistent dim subdim indexing find_order
persistent strict_inputs multimonomial ip
if isempty(dim)
  from labtools import strict_inputs

  from speclab.common.tensor import space_dimension as dim
  from speclab.common.tensor import subspace_dimension as subdim
  from speclab.common.tensor import linear_to_array_indexing as indexing
  from speclab.common.tensor import Npoints_to_poly_order as find_order

  % For default ip, basis:
  from speclab.monomials import multimonomial
  %from speclab.orthopoly.interp import monomial_ip_full as ip
  %from speclab.orthopoly.interp import jacobi_ip as ip
  ip = @(d,k) speye(subdim(d,k));
end

opt = strict_inputs({'tol', 'basis', 'ip'}, {1e-10, [], []}, [], varargin{:});
if isempty(opt.basis)
  bases = cell([size(theta,2) 1]);
  for d = 1:length(bases)
    interval = [min(theta(:,d)) max(theta(:,d))];
    bases{d} = LegendrePolynomialBasis('domain', interval);
  end
  TL = TensorProductBasis(bases{:});
  opt.basis = @(x,n) TL.evaluate(x,n);
end
if isempty(opt.ip)
  opt.ip = ip;
end

[N,d] = size(theta);
l = eye(N);
u = eye(N);
p = speye(N);

% This is just a guess: this vector could be much larger, or much smaller
v = zeros([1000 1]);
v_index = 1;

% Current polynomial degree
k_counter = 0;
k = zeros([N 1]);  % k(q) gives the degree used to eliminate the q'th point

% The current LU row to factor out:
lu_row = 1;
find_rank = true;

while lu_row < N+1
  % We're going to generate the appropriate columns of W -- these are polynomial
  % indices for degree k of the basis.
  current_dim = subdim(d,k_counter);  % The current size of k-vectors
  poly_indices = dim(d,k_counter-1) + (1:current_dim);
  W = p*opt.basis(theta, poly_indices);

  % Row-reduce W according to previous elimination steps
  for q = 1:lu_row-1
    W(q,:) = W(q,:)/l(q,q);
    W(q+1:end,:) = W(q+1:end,:) - l(q+1:end,q)*W(q,:);
  end

  % The mass matrix defining the inner product for this degree
  M = opt.ip(d, k_counter);
  sM = chol(M);

  % Submatrix for qr-like computation:
  temp = full(W(lu_row:N,:)*(sM')).';
  if find_rank
    rnk = rank(temp, opt.tol);
    find_rank = false;
  end
 
  % Find the column with the largest norm, and compute inner products
  norms = sqrt(sum(temp.^2, 1));
  [garbage, next_ind] = max(norms);
  magic_row = temp(:,next_ind).'/norms(next_ind);

  ips = magic_row*temp;
  ips(next_ind) = [];
  ips = [norms(next_ind) ips];

  % Now figure out permutation
  temp_p = speye(length(ips));
  temp_row = temp_p(next_ind,:);
  temp_p(next_ind,:) = [];
  temp_p = [temp_row; temp_p];
  p(lu_row:N,:) = temp_p*p(lu_row:N,:);

  % Permute rows of l:
  l(lu_row:N,1:lu_row-1) = temp_p*l(lu_row:N,1:lu_row-1);

  % Now update l with ip information
  l(lu_row:N,lu_row) = ips;

  % Now find ips with rows above
  u(1:lu_row-1,lu_row) = W(1:lu_row-1,:)*M*magic_row.';

  % Save current information
  v(v_index:v_index+current_dim-1) = magic_row.';

  % Update counters
  v_index = v_index + current_dim;
  k(lu_row) = k_counter;
  lu_row = lu_row + 1;
  if rnk < 2
    k_counter = k_counter + 1;
    find_rank = true;
  else
    rnk = rnk-1;
  end
end

% Chop off parts of unnecessarily allocated vector v
v(v_index:end) = [];
