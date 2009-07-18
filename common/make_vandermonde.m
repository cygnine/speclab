function[v] = make_vandermonde(x,n,functions,varargin)
% [V] = MAKE_VANDERMONDE(X,N,FUNCTIONS,{ARGS})
%
%     The generic Vandermonde matrix constructor for spectral methods. X is the
%     set of nodal points. N is the vector of function indices, and FUNCTIONS is
%     a function handle that takes in arguments X, N, and optionally anything
%     given in the input struct ARGS. Assumes that FUNCTIONS is vectorized in X
%     and N, and that each column of the return value corresponds to each value
%     in the vector N. 

if nargin > 3
  v = functions(x,n,varargin{1});
else
  v = functions(x,n);
end
