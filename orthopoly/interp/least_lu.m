function[l,u,p,v,k_storage] = least_lu(theta, varargin)
% least_lu -- Computes de Boor's LW factorization
%
% [l,u,p,v,k] = least_lu(theta, [tol=1e-10, ip=alpha!,  basis='monomials'])
%
%     Performs the LW factorization as given in [1]. The input theta is
%     an N x D matrix, where there are assumed to be N points, each of which is
%     in D-dimensional space. The result is a factorization such that p*V =
%     l*W, where V is the Vandermonde matrix generated from the basis
%     functions (optional input "basis").  Orthogonalization is performed
%     degree-by-degree using the given inner product (optional input "ip"). The
%     output vector k has length size(theta,1), and for each theta entry, it
%     gives the degree polynomial used to eliminate that row of p*theta.
%
%     The second output u is the N x N upper-triangular matrix that gives inner
%     product information about the elements in W (see [1]). It is necessary to
%     invert u to determine the interpolatory coefficients for multimonomials.
%     Most of the remainder of the matrix W is not necessary; only `diagonal'
%     elements of it are stored, and are given in the vector v.
%
%     The optional input tol is the tolerance for the norm of a potential pivot
%     below which it is assumed to vanish. 
%
%     The input "basis" should be a function handle with the following syntax:
%
%       vandermonde_matrix = basis(theta, 0:(M-1))
%
%     which generates a size(theta,1) x M matrix where column m is the
%     evaluation of basis function m at the multidimensional locations theta.
%
%     The input "ip" should be a function handle with the following syntax:
%
%       M = ip(d, k)
%
%     which generates a square matrix M defining the inner product between two
%     vectors representing the coefficients of the basis functions in d-variate
%     polynomial space of total degree equal to k. The size of M should be
%     dictated by the size of this space.
%
%     The default values of `basis' and `ip' generate the "least polynomial" of
%     de Boor and Ron.
% 
% [1]:  Computational Aspects of Polynomial Interpolation in Several Variables,
%       Carl de Boor and Amos Ron, Mathematics of Computation, 1992.

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
  from speclab.orthopoly.interp import monomial_ip as ip
end

opt = strict_inputs({'tol', 'basis', 'ip'}, {1e-10, [], []}, [], varargin{:});
if isempty(opt.basis)
  opt.basis = @(x,n) multimonomial(x,n,'dim',size(theta,2));
end
if isempty(opt.ip)
  opt.ip = ip;
end

[N,D] = size(theta);
l = eye(N);
u = eye(N);
p = speye(N);
stemp = spalloc(1, N, 2);

% This is just a guess: this vector could be much larger, or much smaller
v = zeros([1000 1]);
v_index = 1;

% Current polynomial degree
k_counter = 0;
k_storage = zeros([N 1]);  % k_storage(q) gives the degree used to eliminate the q'th point

% The current LU row to factor out:
lu_row = 1;

while lu_row < N+1
  current_dim = subdim(D, k_counter); % The current size of k-vectors
  cols = dim(D,k_counter-1)+1;
  cols = cols:(cols+current_dim-1); % A vector to pick out degree-k columns from W

  % First generate appropriate columns of Vandermonde matrix
  W = p*opt.basis(theta, cols-1);

  % Use diary of operations to bring to upper triangular form
  for q = 1:(lu_row-1)
    W(q+1:end,:) = W(q+1:end,:) - l(q+1:end,q)*W(q,:);
  end

  % The mass matrix defining the inner product for this degree
  M = opt.ip(D, k_counter);

  % Force one elimination: in exact arithmetic, one is always needed
  first_pass = true;
  another_pass = true;

  % Compute degree-k norms for pivoting
  norms = sum(abs(W(lu_row:end, :)).^2*M, 2);

  while another_pass

    % Just a clause to save some minimal computation if we're at the last entry
    if lu_row == N
      % This is the case when we're at the last row, and it has a nontrivial
      % norm. So just calculate appropriate stuff and break and return.
      u(1:lu_row, lu_row) = W(1:lu_row,:)*M*W(lu_row,:).';
      v(v_index:v_index+current_dim-1) = W(lu_row,:);
      v_index = v_index + current_dim;
      k_storage(lu_row) = k_counter;
      lu_row = lu_row+1;
      break;
    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
    %%%%%%%%%%    PIVOTING    %%%%%%%%%% 
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
    
    [garbage, pivot] = max(norms);
    pivot = pivot + lu_row - 1;

    % When pivoting, must rearrange W, l, and p
    vtemp = W(pivot,:);
    W(pivot,:) = W(lu_row,:);
    W(lu_row,:) = vtemp;

    vtemp = l(pivot,1:lu_row-1);
    l(pivot,1:lu_row-1) = l(lu_row,1:lu_row-1);
    l(lu_row,1:lu_row-1) = vtemp;

    stemp(:,:) = p(lu_row,:);
    p(lu_row,:) = p(pivot,:);
    p(pivot,:) = stemp;

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
    %%%%%%%%%%  END PIVOTING  %%%%%%%%%% 
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 

    % Calculate the inner products and store in the upper-triangular matrix u:
    u(1:lu_row, lu_row) = W(1:lu_row,:)*M*W(lu_row,:).';

    % Gauss eliminate:
    factors = W(lu_row,:)*M*W(lu_row+1:end,:).';
    factors = factors/u(lu_row,lu_row);
    l(lu_row+1:end,lu_row) = factors.';
    W(lu_row+1:end,:) = W(lu_row+1:end,:) - factors.'*W(lu_row,:);

    % Store `diagonal' elements of W
    v(v_index:v_index+current_dim-1) = W(lu_row,:);
    v_index = v_index + current_dim;

    % Update values
    k_storage(lu_row) = k_counter;
    lu_row = lu_row + 1;

    % See if any rows below have nontrivial norm
    norms = sum(abs(W(lu_row:end, :)).^2*M, 2);
    first_pass = false;
    another_pass = first_pass | (any(norms>opt.tol) & (lu_row < N+1));
  end

  k_counter = k_counter + 1;
end

% Chop off unnecessarily allocated v
v(v_index:end) = [];
