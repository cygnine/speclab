function[y] = stiffness_operator(x,varargin)
% [y] = stiffness_operator(x,{alpha=-1/2,shift=0,scale=1,normalization='normal'})
%
%     Uses an O(N) algorithm to compute modal coefficients y for the derivative of
%     a modal expansion given the modal coefficients x of the function. The
%     optional parameters help specifiy the particular class of Laguerre
%     polynomials.
%
%     This method can be seen as an implementation of one of two equivalent
%     algorithms:
%       - Determining and solving the three-term backward recurrence for doing Laguerre
%         derivatives.
%       - Using coefficients to find the derivative expansion in terms of
%         polynomials of class (alpha+1). Then use the sparse integer-separation
%         connection to invert a tridiagonal upper-diagonal system. 

persistent defaults matinv derivative connection
if isempty(defaults)
  from speclab.orthopoly.laguerre import defaults
  from labtools.linalg import triu_sparse_invert as matinv
  from speclab.orthopoly.laguerre import derivative_connection as derivative
  from speclab.orthopoly.laguerre import integer_separation_connection_matrix as connection
end

opt = defaults(varargin{:});

% Force column vector
x = x(:);
N = length(x);

% Get derivative coefficients
zetas = derivative(0:(N-1),opt.alpha,opt);
zetas = zetas*opt.scale;  % WTF!?!!?!

% Get connection matrix
C = connection(N,opt.alpha,1,opt);

% new coefficients in (alpha+1)
y = x(2:end).*zetas(2:end);
y = [y;0];

y = matinv(C,y,'bandwidth',3);
