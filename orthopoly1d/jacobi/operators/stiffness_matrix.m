function[S] = stiffness_matrix(N,varargin)
% [S] = stiffness_matrix(N,{alpha=-1/2,beta=-1/2,shift=0,scale=1,normalization='normal'})
%
%     Returns the N x N stiffness matrix for the Jacobi polynomials of class
%     (alpha, beta). Since the matrix has at most ~(N^2)/2 nonzero elements,
%     it's returned as a sparse matrix.
%
%     This is the direct method for stiffness_operator.

global packages;
jac = packages.speclab.orthopoly1d.jacobi;
opt = jac.defaults(varargin{:});

% Force column vector
S = zeros(N);

% Get derivative coefficients
zetas = jac.coefficients.derivative(0:(N-1),opt.alpha,opt.beta,opt);
zetas = zetas/opt.scale;  % WTF!?!!?!

% Get connection matrix
C = jac.connection.integer_separation_connection_matrix(N,opt.alpha,opt.beta,1,1,opt);

% new coefficients in (alpha+1,beta+1):
S = spdiags(zetas,1,N,N);
S = inv(C)*S;