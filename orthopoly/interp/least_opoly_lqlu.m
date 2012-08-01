function[varargout] = least_opoly_lqlu(theta, varargin)
%function[l,u,p,v,k,H] = least_opoly_lqlu(theta, varargin)
% tensor_lu -- Computes de Boor's LU factorization with LQ factorizations
%
% [L,U,P,H] = least_opoly_lqlu(theta, [tol=1e-10, ip=eye(),  basis='legendre'])
% [L,U,P,H,v,k] = least_opoly_lqlu(theta, [tol=1e-10, ip=eye(),  basis='legendre'])
%
%     Performs the LUH factorization for the least orthogonal interpolant as
%     given in [1]. The input theta is an N x D matrix, where there are assumed
%     to be N points, each of which is in D-dimensional space. Given an input
%     polynomial basis, this makes use of a Vandermonde-like matrix V with
%     entries 
%
%               V(n,m) = basis_m(theta(n,:)) 
%
%     where basis_m means the m'th basis coefficients. The least orthogonal
%     space is formed by a tensor-product Legendre basis, where the univariate
%     basis in dimension d is formed by the smallest enclosing interval for
%     theta(:,d).
%
%     OUTPUTS: 
%     The output of this function is a matrix factorization of V:
%
%               P*V ~= L*U*H,
%
%     where P is an N x N permutation matrix, L is an N x N lower-triangular
%     matrix, U is an N x N upper triangular matrix, and H is an N x R matrix,
%     where R is the number of columns of V, which is not known a priori. The
%     matrices L, U, and P are the same as in a row-pivoting L U
%     factorization.  The matrix factorization is not exact because irrelevant
%     entries in H are discarded resulting in a block-diagonal matrix. 
%
%     The matrix H is effectively a basis defining the least orthgonal space.
%     If the least orthogonal space is a full-degree polynomial space, then the
%     matrix factorization is an equality. 
%
%     The output v is simply a linear vector containing entries of H, and the
%     output k is a vector of length N, with entry n indicating the polynomial
%     degree associated with input point n. k has non-decreasing entries.
%
%     The least orthogonal pseudoinverse action of V, which maps data f to
%     input basis coefficients c is given by
%
%               c = H'*inv(L*U)*P*f
%
%     OPTIONAL INPUTS:
%     ip is a function handle with syntax ip(k), which returns the Gram matrix
%     for degree k. (For orthonormal polynomials, this is the identity matrix.)
%
%     tol is a scalar tolerance for 2-norms below which entries are treated as
%     0.
%
%     basis can be given to change the least orthogonal space (by providing
%     e.g. tensor-product Hermite polynomials). It is either a
%     TensorProductBasis instance, or a function handle supporting the calling
%     syntax basis(theta,m), which for N x D matrix theta, and length-M vector
%     m, produces an N x M Vandermonde-like matrix V. The (j,k) entry of V
%     should be the k'th basis function evaluated at theta(j,:). (The basis
%     should be total-degree ordered and 1-indexed, k >= 1.)

persistent dim subdim spdiag
persistent parser ip
if isempty(dim)
  from labtools import input_parser spdiag

  from speclab.common.tensor import polynomial_space_dimension as dim
  from speclab.common.tensor import polynomial_subspace_dimension as subdim

  [opt, parser] = input_parser({'tol', 'basis', 'ip'}, ...
                               {1e-10, [], []}, ...
                               [], ...
                               varargin{:});

  % For default ip, basis:
  ip = @(d,k) speye(subdim(d,k));
else
  parser.parse(varargin{:});
  opt = parser.Results;
end

if isempty(opt.basis)
  bases = cell([size(theta,2) 1]);
  for d = 1:length(bases)
    interval = [min(theta(:,d)) max(theta(:,d))];
    bases{d} = LegendrePolynomialBasis('domain', interval);
  end
  %TL = TensorProductBasis(bases, 'indexing', 1);
  TL = TensorProductBasis(bases);
  opt.basis = @(x,n) TL.evaluate(x,n);
end
if isempty(opt.ip)
  opt.ip = ip;
end

[N,d] = size(theta);
l = eye(N);
u = eye(N);
p = speye(N);
%stemp = spalloc(1, N, 2);

% This is just a guess: this vector could be much larger, or much smaller
v = zeros([1000 1]);
v_index = 1;

% Current polynomial degree
k_counter = 0;
k = zeros([N 1]);  % k(q) gives the degree used to eliminate the q'th point

% The current LU row to factor out:
lu_row = 1;

% Current degree is k_counter, and we iterate on this
while lu_row < N+1
  % We're going to generate the appropriate columns of W -- these are polynomial
  % indices for degree k of the basis.
  current_dim = subdim(d,k_counter);  % The current size of k-vectors

  if isa(opt.basis, 'TensorProductBasis')
    poly_indices = opt.basis.range(dim(d,k_counter-1) + current_dim);
    poly_indices = poly_indices((end-current_dim+1):end, :);
    %%%%%%%%% PRECONDITIONING %%%%%%%%%%%
    W = p*full(spdiag(sqrt(opt.basis.weight(theta)))*opt.basis(theta, poly_indices));
  else
    % Assume 1-based indexing
    poly_indices = dim(d,k_counter-1) + (1:current_dim);
    W = p*opt.basis(theta, poly_indices);
  end

  % Row-reduce W according to previous elimination steps
  for q = 1:lu_row-1
    W(q,:) = W(q,:)/l(q,q);
    W(q+1:end,:) = W(q+1:end,:) - l(q+1:end,q)*W(q,:);
  end

  % The mass matrix defining the inner product for this degree
  M = opt.ip(d, k_counter);
  sM = chol(M);

  % Determine submatrix indices required for qr factorization
  %[q,r,e] = qr((W(lu_row:N,:)*(sM')).');
  [q,r,evec] = qr(full(W(lu_row:N,:)*(sM')).', 0);
  rnk = rank(r, opt.tol);

  NN = length(lu_row:N);
  e = spalloc(NN,NN,NN);
  for qq = 1:NN
    e(evec(qq),qq) = 1;
  end

  % Now first we must permute the rows by e
  p(lu_row:N,:) = sparse(e')*p(lu_row:N,:);
  % And correct by permuting l as well:
  l(lu_row:N,1:lu_row-1) = e'*l(lu_row:N,1:lu_row-1);
  
  % The matrix r gives us inner product information for all rows below these in
  % W
  l(lu_row:N,lu_row:lu_row+rnk-1) = r(1:rnk,:).';

  % Now we must find inner products of all the other rows above these in W
  u(1:lu_row-1,lu_row:lu_row+rnk-1) = W(1:lu_row-1,:)*M*q(:,1:rnk);

  % The matrix q must be saved in order to characterize basis
  v(v_index:v_index+(current_dim*rnk)-1) = reshape(q(:,1:rnk), [rnk*current_dim 1]);
  v_index = v_index+(current_dim*rnk);

  % Update degree markers, and node and degree count
  k(lu_row:(lu_row+rnk-1)) = k_counter;
  lu_row = lu_row + rnk;
  k_counter = k_counter + 1;
end

% Chop off parts of unnecessarily allocated vector v
v(v_index:end) = [];

% Make matrix H:
H = spalloc(N,k_counter-1,length(v));
v_counter = 1;
for q = 1:N;
  current_size = subdim(d, k(q));
  v_inds = v_counter:(v_counter+current_size-1);
  
  previous_dimension = dim(d, k(q)-1);
  H_cols = (previous_dimension+1):(previous_dimension+current_size);

  H(q,H_cols) = v(v_inds);

  v_counter = v_counter + current_size;
end

switch nargout
case 4
  varargout{1} = l;
  varargout{2} = u;
  varargout{3} = p;
  varargout{4} = H;
case 6
  varargout{1} = l;
  varargout{2} = u;
  varargout{3} = p;
  varargout{4} = H;
  varargout{5} = v;
  varargout{6} = k;
end
