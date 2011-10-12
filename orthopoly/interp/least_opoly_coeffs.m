function[c] = least_opoly_coeffs(theta, f, varargin)
% least_opoly_coeffs -- Computes basis coefficients for interpolation
%
% c = least_opoly_coeffs(l,u,p,v,k,dim,f, varargin)
%
%     Basically stitches together least_opoly_lqlu and least_opoly_lucoffs --
%     takes in data points (theta, f) and returns coefficients c from the
%     specified basis. The optional inputs are those that can be given to
%     least_opoly_lqlu and least_opoly_lucoeffs.
%
%     The interpolant may be evaluated as 
%
%       p(x) = sum_n c(n) * p(x,n),
%
%     where p(x,n) is a function that evaluates the n'th basis function at the
%     location x. (e.g. speclab.monomials.multimonomial is one such function). 
%
%     Note that this function is vectorized in the columns of f, meaning that if
%     f is a collection of column vectors, this function returns a collection of
%     column vectors containing the basis representation.

persistent opoly_lu opoly_coeffs
if isempty(opoly_lu)
  from speclab.orthopoly.interp import least_opoly_lqlu as opoly_lu
  from speclab.orthopoly.interp import least_opoly_lucoeffs as opoly_coeffs
end

[l,u,p,v,k] = opoly_lu(theta, varargin{:});
c = opoly_coeffs(l,u,p,v,k,size(theta,2),f, varargin{:});
