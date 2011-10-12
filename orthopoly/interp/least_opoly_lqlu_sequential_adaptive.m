function[l,u,p,v,k,inds] = least_opoly_lqlu_sequential_adaptive(theta, f, N, varargin)
% least_opoly_lqlu_sequential_adaptive -- Computes de Boor's LU factorization with LQ factorizations
%
% [l,u,v,p,k] = least_opoly_lqlu_sequential_adaptive(theta, f, N, [tol=1e-10, ip=eye(),  basis='legendre', initial_inds = 1])
%
%     Using the sequential method, this routine chooses N points using the
%     (theta,f) data that minimizes the change in the coefficients. 
%
%     We always start out with an initial collection of nodes. The matrix
%     initial_inds determines this.

persistent dim subdim indexing find_order spdiag
persistent strict_inputs multimonomial ip invu
if isempty(dim)
  from labtools import strict_inputs
  from labtools import spdiag

  from speclab.common.tensor import space_dimension as dim
  from speclab.common.tensor import subspace_dimension as subdim
  from speclab.common.tensor import linear_to_array_indexing as indexing
  from speclab.common.tensor import Npoints_to_poly_order as find_order
  from labtools.linalg import triu_back_substitute as invu

  % For default ip, basis:
  from speclab.monomials import multimonomial
  ip = @(d,k) speye(subdim(d,k));
end

opt = strict_inputs({'tol', 'basis', 'ip', 'initial_inds'}, {1e-10, [], [], 1}, [], varargin{:});
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

[Nt,d] = size(theta);

% Can't adapt more points than we're given 
if N > Nt
  N = Nt;
end

%l = eye(N);
u = eye(N);
p = speye(Nt);

% This stores inner product information for all nodes -- we'll need this in
% order to update coefficients once we choose points. It's faster than
% recomputing them.
ltemp = zeros([Nt N]);

% This is just a guess: this vector could be much larger, or much smaller
v = zeros([1000 1]);
v_index = 1;

% Current polynomial degree
k = zeros([N 1]);  % k(q) gives the degree used to eliminate the q'th point

% The current LU row to factor out:
find_rank = true;

% Arggggghhhh...find permutation matrix to shift initial_inds to the top.
Ninit = length(opt.initial_inds);
%temp = p(opt.initial_inds,:);
%p(opt.initial_inds,:) = [];
%p = [temp; p];

lu_row = 1;
k_counter = 0;
v_index = 1;

cprev = 0;

while lu_row < N+1
  % We're going to generate the appropriate columns of W -- these are polynomial
  % indices for degree k of the basis.
  current_dim = subdim(d,k_counter);  % The current size of k-vectors
  poly_indices = dim(d,k_counter-1) + (1:current_dim);
  W = p*opt.basis(theta, poly_indices);

  % Row-reduce W according to previous elimination steps
  for q = 1:lu_row-1
    W(q,:) = W(q,:)/ltemp(q,q);
    W(q+1:end,:) = W(q+1:end,:) - ltemp(q+1:end,q)*W(q,:);
  end

  % The mass matrix defining the inner product for this degree
  M = opt.ip(d, k_counter);
  sM = chol(M);

  % Submatrix for qr-like computation:
  temp = full(W(lu_row:Nt,:)*(sM')).';
  if find_rank
    rnk = rank(temp, opt.tol);
    if rnk < 1; % Oops, we actually can't use this degree
      k_counter = k_counter + 1;
      continue;
    else
      find_rank = false;
    end
  end

  if lu_row <= Ninit
    % Then it's easy:
    next_ind = opt.initial_inds(lu_row) - lu_row + 1;
  else % Need to do some work
    % Now for each points lu_row ..... Nt, we want to (a) compute the next
    % forward-substitution coefficient, and (b) choose the one with the smallest
    % forward-substitution coefficient.

    % First: what would the norms be if lu_row, ..., Nt were the next point?
    norms = sqrt(sum(temp.^2,1).');

    % But we must exclude points that force a new degree
    bad_inds = norms <= opt.tol;

    %candidate_coeffs = ones([size(temp,2) 1])*f(lu_row) - ltemp(lu_row:Nt,1:lu_row-1)*f(1:lu_row-1);
    candidate_coeffs = f(lu_row:Nt) - ltemp(lu_row:Nt,1:lu_row-1)*f(1:lu_row-1);
    candidate_coeffs = candidate_coeffs./norms;

    % Cheap hack fix
    candidate_coeffs(bad_inds) = Inf;

    % Compute 'upper' inner products
    utemp = W(1:lu_row-1,:)*M*temp*spdiag(1./norms);

    % Now back substitute 
    ftemp = repmat(f(1:lu_row-1), [1 size(temp,2)]) - utemp*spdiag(candidate_coeffs);
    %ftemp = inv(u(1:lu_row-1,1:lu_row-1))*ftemp;
    %ftemp = uinv*ftemp;
    ftemp = invu(u(1:lu_row-1,1:lu_row-1), ftemp);

    candidate_norms = ftemp - repmat(cprev, [1 size(temp,2)]);
    %candidate_norms = ftemp;
    candidate_norms = sum(candidate_norms.^2,1).' + candidate_coeffs.^2;
    %candidate_norms = sum(abs(candidate_norms),1).' + abs(candidate_coeffs);

    % Now instead we're going to minimize the coefficient difference.
    %[garbage, next_ind] = min(abs(candidate_coeffs));
    [garbage, next_ind] = min(abs(candidate_norms));
    %fprintf('%3d  %1.3e\n', lu_row, garbage);
  end

  % Now compute inner products, etc.
  norms = norm(temp(:,next_ind));
  magic_row = temp(:,next_ind).'/norms;

  ips = magic_row*temp;
  ips(next_ind) = [];
  ips = [norms ips];

  %if abs(ips(1)) < 1e-9
  %  asdf
  %end

  % Now figure out permutation
  temp_p = speye(length(ips));
  temp_row = temp_p(next_ind,:);
  temp_p(next_ind,:) = [];
  temp_p = [temp_row; temp_p];
  p(lu_row:Nt,:) = temp_p*p(lu_row:Nt,:);

  f(lu_row:Nt) = temp_p*f(lu_row:Nt);

  % Permute rows of l:
  ltemp(lu_row:Nt,1:lu_row-1) = temp_p*ltemp(lu_row:Nt,1:lu_row-1);

  % Now update l with ip information
  ltemp(lu_row:Nt,lu_row) = ips;

  % Now find ips with rows above
  u(1:lu_row-1,lu_row) = W(1:lu_row-1,:)*M*magic_row.';

  % Save current information
  v(v_index:v_index+current_dim-1) = magic_row.';

  % Update rhs vector:
  f(lu_row) = f(lu_row) - ltemp(lu_row,1:lu_row-1)*f(1:lu_row-1);
  f(lu_row) = f(lu_row)/ltemp(lu_row,lu_row);
  %uinv = inv(u(1:lu_row,1:lu_row));
  %cprev = uinv*f(1:lu_row);
  cprev = invu(u(1:lu_row,1:lu_row), f(1:lu_row));

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

l = ltemp(1:N,1:N);
inds = (1:Nt).';
inds = p(1:N,:)*inds;
p = speye(N);

%f(1:N);
%cprev;

% Chop off parts of unnecessarily allocated vector v
v(v_index:end) = [];
