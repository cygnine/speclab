function[y] = stiffness_operator(x,varargin)
% [Y] = STIFFNESS_OPERATOR(X,{ALPHA=-1/2,BETA=-1/2,SHIFT=0,SCALE=1,NORMALIZATION='NORMAL'})
%
%     Uses an O(N) algorithm to compute modal coefficients Y for the derivative of
%     a modal expansion given the modal coefficients X of the function. The
%     optional parameters help specifiy the particular class of Jacobi
%     polynomials.
%
%     This method can be seen as an implementation of one of two equivalent
%     algorithms:
%       - Determining and solving the three-term backward recurrence for doing Jacobi
%         derivatives.
%       - Using coefficients to find the derivative expansion in terms of
%         polynomials of class (ALPHA+1,BETA+1). Then use the sparse integer-separation
%         connection to invert a tridiagonal upper-diagonal system. 

global handles;
jac = handles.speclab.OrthogonalPolynomial1D.jacobi;
opt = jac.defaults(varargin{:});

% Force column vector
x = x(:);
N = length(x);

% Get derivative coefficients
zetas = jac.coefficients.derivative(0:(N-1),opt);

% Get connection matrix
C = jac.connection.integer_separation_connection_matrix(N,opt.alpha,opt.beta,1,1,opt);

% new coefficients in (alpha+1,beta+1):
y = x(2:end)./zetas(2:end);
y = [y;0];

%%% GET SPARSE INVERSION LA UP DUMMY %%%%
error('Not implemented yet');
