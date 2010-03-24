function[l,W,p,k_storage,u] = opoly_li_lu(theta, varargin)
% opoly_li_lu -- Computes de Boor's LW factorization for orthoogonal polynomials
%
% [l,W,p,k,u] = opoly_li_lu(theta, [tol=1e-10, opoly_class='jacobi', 
%                                   opoly_params={'alpha', 0, 'beta', 0}])
%
%     Performs the LW factorization as given in [1] using (1) an orthogonal
%     polynomial basis and (2) the inner product defined by the L^2
%     orthogonality relation defined by the polynomial class. The input theta is
%     an N x D matrix, where there are assumed to be N points, each of which is
%     in D-dimensional space. The result is a factorization such that p*V = l*W,
%     where V is the orthogonal polynomial Vandermonde matrix generated in the
%     code. The output vector k has length size(theta,1), and for each theta
%     entry, it gives the degree polynomial used to eliminate that row of
%     p*theta.
%
%     The final output u is the N x N upper-triangular matrix that gives inner
%     product information about the elements in W (see [1]). It is necessary to
%     invert u to determine the interpolatory coefficients for orthogonal
%     polynomials.
%
%     The optional input tol is the tolerance for the norm of a potential pivot
%     below which it is assumed to vanish. The norm is measured using the
%     associated inner product.
% 
% [1]:  Computational Aspects of Polynomial Interpolation in Several Variables,
%       Carl de Boor and Amos Ron, Mathematics of Computation, 1992.

persistent dim subdim indexing find_order
persistent eval_jacobi_poly eval_hermite_poly
persistent ip strict_inputs
if isempty(dim)
  from labtools import strict_inputs

  from speclab.common.tensor import space_dimension as dim
  from speclab.common.tensor import subspace_dimension as subdim
  from speclab.common.tensor import linear_to_array_indexing as indexing
  from speclab.common.tensor import Npoints_to_poly_order as find_order
  from speclab.orthopoly1d.jacobi.eval import eval_jacobi_poly
  from speclab.orthopoly1d.hermite.eval import eval_hermite_poly

  ip = @(d,k) speye(subdim(d,k)); % L^2 orthonormal polynomials are nice
end

opt = strict_inputs({'tol', 'opoly_class', 'opoly_params'}, ...
                     {1e-10, 'jacobi', struct('alpha', 0, 'beta', 0)}, [], varargin{:});

switch lower(opt.opoly_class)
case {'jacobi', 'jac'}
  poly_eval = eval_jacobi_poly;
case {'hermite', 'herm'}
  poly_eval = eval_hermite_poly;
end

[N,D] = size(theta);
l = eye(N);
u = eye(N);
p = speye(N);
stemp = spalloc(1, N, 2);

% Current polynomial degree
k_counter = 0;
k_storage = zeros([N 1]);  % k_storage(q) gives the degree used to eliminate the q'th point

% The current LU row to factor out:
lu_row = 1;

% Temporary variables to compute the initial Vandermonde matrix:
highest_degree = find_order(D, N);
poly_indices = 0:(dim(D, highest_degree)-1);

% Generate initial Vandermonde matrix:
opt.opoly_params.dim = D;
W = poly_eval(theta, poly_indices, opt.opoly_params);
[Wr,Wc] = size(W);
% If during the algorithm we need more columns than Wc, they will be generated
% at that point.

while lu_row < N+1
  % The mass matrix defining the inner product for this degree
  M = ip(D, k_counter);

  current_dim = subdim(D, k_counter); % The current size of k-vectors
  cols = dim(D,k_counter-1)+1;
  cols = cols:(cols+current_dim-1); % A vector to pick out degree-k columns from W

  % Force one elimination: in exact arithmetic, one is always needed
  first_pass = true;
  another_pass = true;

  % Compute degree-k norms for pivoting
  norms = sum(abs(W(lu_row:end, cols)).^2*M, 2);

  while another_pass

    if lu_row == N
      % This is the case when we're at the last row, and it has a nontrivial
      % norm. So just calculate appropriate stuff for u and break and return.

      % Calculate the norms for the upper-triangular matrix u:
      u(1:lu_row, lu_row) = W(1:lu_row,cols)*M*W(lu_row,cols).';

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

    % Gauss eliminate:
    factors = W(lu_row,cols)*M*W(lu_row+1:end,cols).';
    factors = factors/(W(lu_row,cols)*M*W(lu_row,cols).');
    l(lu_row+1:end,lu_row) = factors.';
    W(lu_row+1:end,cols(1):end) = W(lu_row+1:end,cols(1):end) - factors.'*W(lu_row,cols(1):end);

    % Calculate the norms for the upper-triangular matrix u:
    u(1:lu_row, lu_row) = W(1:lu_row,cols)*M*W(lu_row,cols).';

    k_storage(lu_row) = k_counter;
    lu_row = lu_row + 1;

    norms = sum(abs(W(lu_row:end, cols)).^2*M, 2);

    first_pass = false;

    another_pass = first_pass | (any(norms>opt.tol) & (lu_row < N+1));
  end

  k_counter = k_counter + 1;

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %%%%%%%%%%%% EXTEND W %%%%%%%%%%%%%%
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  if (highest_degree < k_counter) & (lu_row < N+1)
    % 1.) Generate a temporary Vandermonde matrix for the extra degree. 
    % 2.) Permute as given by p
    % 3.) Gauss-eliminate as determined by l
    % 4.) Append to current W

    poly_indices = dim(D,k_counter-1);
    poly_indices = poly_indices:(dim(D,k_counter)-1);
    Wtemp = poly_eval(theta, poly_indices, opt.opoly_params);

    Wtemp = p*Wtemp;

    for q = 1:(lu_row-1)
      Wtemp(q+1:end,:) = Wtemp(q+1:end,:) - l(q+1:end,q)*Wtemp(q,:);
    end

    W = [W Wtemp];

    highest_degree = highest_degree + 1;
  end
end
