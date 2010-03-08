function[fz] = monomial_interp(theta, f, z, varargin)
% monomial_interp -- Computes multidimensional interpolation via de Boor's algorithm
%
% fz = monomial_interp(theta, f, z, [tol=1e-10])
%
%     Using de Boor's algorithm [1] this function determines the `least
%     interpolant' using the location-data pair (theta, f). The dimension of the
%     space is automatically set to size(theta, 2). f and theta must have the
%     same number of rows. z is another collection of points in space at which
%     the least interpolant is evaluated, which is returned in fz.
%
%     This function is hardly the most effective way to compute the interpolant
%     -- it is just a patching together of different routines.
%
%     Note that this function is vectorized in the columns of f -- meaning that
%     for each column of f, an interpolant is constructed and evaluated at the
%     locations z. Therefore fz has as many rows as z and as many columns as f.

persistent monomial_lu monomial_coeffs multimonomial
persistent strict_inputs
if isempty(strict_inputs)
  from labtools import strict_inputs
  
  from speclab.orthopoly1d.interp import monomial_lu monomial_coeffs
  from speclab.monomials import multimonomial
end

opt = strict_inputs({'tol'}, {1e-10}, [], varargin{:});

[l,W,p,k,u] = monomial_lu(theta, 'tol', opt.tol);
dim = size(theta, 2);

c = monomial_coeffs(l,W,p,k,u,dim,f);

V = multimonomial(z, 0:(length(c)-1), 'dim', dim);

fz = V*c;
