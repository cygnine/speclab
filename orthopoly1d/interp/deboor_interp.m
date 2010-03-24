function[fz] = deboor_interp(theta, f, z, varargin)
% deboor_interp -- Computes multidimensional interpolation via de Boor's algorithm
%
% fz = deboor_interp(theta, f, z, [tol=1e-10, ip=alpha!, basis=multimonomial])
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

persistent deboor_lu deboor_coeffs multimonomial monomial_ip
persistent strict_inputs
if isempty(strict_inputs)
  from labtools import strict_inputs
  
  from speclab.orthopoly1d.interp import deboor_lu deboor_coeffs
  from speclab.orthopoly1d.interp import monomial_ip
  from speclab.monomials import multimonomial
end

opt = strict_inputs({'tol','ip','basis'}, {1e-10, monomial_ip, []}, [], varargin{:});

dim = size(theta, 2);
if isempty(opt.basis)
  opt.basis = @(x,n) multimonomial(x,n,'dim',dim);
end

[l,W,p,k,u] = deboor_lu(theta, opt);

c = deboor_coeffs(l,W,p,k,u,dim,f,'ip',opt.ip,'tol',opt.tol);

V = opt.basis(z, 0:(length(c)-1));

fz = V*c;
