function[l,u,p,v,k] = least_opoly_lqlu(theta, varargin)
% tensor_lu -- Computes de Boor's LU factorization with LQ factorizations
%
% [l,u,v,p,k] = least_opoly_lqlu(theta, [tol=1e-10, ip=eye(),  basis='legendre'])
%
%     Performs the LU factorization as given in [1]. The input theta  is
%     an N x D matrix, where there are assumed to be N points, each of which is
%     in D-dimensional space. The result is a factorization such that p*V =
%     l*W, where V is the Vandermonde matrix generated from the basis functions
%     (optional input "basis").  Orthogonalization is performed degree-by-degree
%     using the given inner product (optional input "ip") and a QR
%     factorization. The output vector k has length size(theta,1), and for each
%     theta entry, it gives the degree polynomial used to eliminate that row of
%     p*theta.
%
%     The final output u is the N x N upper-triangular matrix that gives inner
%     product information about the elements in W (see [1]). It is necessary to
%     invert u to determine the interpolatory coefficients for multimonomials.
%
%     The optional input tol is the tolerance for the norm of a potential pivot
%     below which it is assumed to vanish. 
%
%     The input "basis" should be a function handle with the following syntax:
%
%       vandermonde_matrix = basis(theta, 1:M)
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
%     The default values of "basis" and "ip" generate the "least" polynomial of
%     de Boor and Ron.
% 
% [1]:  Computational Aspects of Polynomial Interpolation in Several Variables,
%       Carl de Boor and Amos Ron, Mathematics of Computation, 1992.

persistent dim subdim indexing find_order
persistent parser multimonomial ip
if isempty(dim)
  from labtools import input_parser

  from speclab.common.tensor import space_dimension as dim
  from speclab.common.tensor import subspace_dimension as subdim
  from speclab.common.tensor import linear_to_array_indexing as indexing
  from speclab.common.tensor import Npoints_to_poly_order as find_order

  [opt, parser] = input_parser({'tol', 'basis', 'ip'}, ...
                               {1e-10, [], []}, ...
                               [], ...
                               varargin{:});

  % For default ip, basis:
  from speclab.monomials import multimonomial
  ip = @(d,k) speye(subdim(d,k));
else
  parser.parse(varargin{:});
  opt = parser.Results;
end

if isempty(opt.basis)
  %opt.basis = @(x,n) multimonomial(x,n,'dim',size(theta,2));
  bases = cell([size(theta,2) 1]);
  for d = 1:length(bases)
    interval = [min(theta(:,d)) max(theta(:,d))];
    bases{d} = LegendrePolynomialBasis('domain', interval);
    %bases{d} = LegendrePolynomialBasis();
  end
  TL = TensorProductBasis(bases, 'indexing', 1);
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

  k(lu_row:(lu_row+rnk-1)) = k_counter;

  lu_row = lu_row + rnk;
  k_counter = k_counter + 1;
end

% Chop off parts of unnecessarily allocated vector v
v(v_index:end) = [];
