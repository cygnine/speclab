function[l,u,p] = lu_factorization(theta, varargin)
% lu_factorization -- Determines and factorizes the Vandermonde matrix
%
% [l,u,p] = lu_factorization(theta, [scale=1, shift=0, opoly='jacobi', alpha=0,beta=0)
%
%     Using the D-dimensional data theta (which is N x D), this function returns
%     an LU factorization of the N x N Vandermonde matrix for interpolation at
%     the locations theta. 
%
%     The locations theta are assumed to lie inside the D-dimensional hypercube
%     [-1,1]^D. The optional inputs scale and shift may be used to augment this
%     hypercube. 
%
%     The basis functions for the Vandermonde matrix are tensor-product
%     orthogonal polynomial functions of the specified type (opoly, alpha,
%     beta). At this point, only Jacobi polynomials are implemented. The degree
%     is increased until a suitable (i.e. non-numerically-singular) matrix is
%     found.

persistent strict_inputs npoly_degree
if isempty(strict_inputs)
  from labtools import strict_inputs
  from speclab.common.tensor import npoly_degree
end

opt = strict_inputs({'alpha', 'beta', 'shift', 'scale'}, {0,0,0,1}, [], varargin{:});

singularity_tol = 1e-13;

[N,D] = size(theta);
opt.dim = D;

% First generate Vandermonde matrix (sequence of evaluations):
V = eval_jpoly(theta, 0:ceil(3/2*N), opt);

current_degree = 0;  % polynomial degree counter
n_points = 0;  % Number of points that have been factorized

while n_points < N
  qr_size = npoly_degree(current_degree, D);  % size of the block to factorize
  qr_size = min([qr_size N-n_points]);

  % Now for each block, determine the next qr_size theta values that
  % qr-factorize well.
  l = 
