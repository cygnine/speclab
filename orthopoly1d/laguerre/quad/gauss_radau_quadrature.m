function[x,w] = gauss_radau_quadrature(N,varargin)
% gauss_radau_quadrature -- Gauss-Radau quadrature rule for Laguerre polynomials
%
% [x,w] = gauss_radau_quadrature(N,{alpha=0,shift=0,scale=1,r=0})
%
%     Returns the N-point Gauss-Radau quadrature rule for the Laguerre polynomials.
%     The weight function is given by speclab.orthopoly1d.laguerre.weights.weight.
%     The Radau point is located at x=r.

global packages;
opoly = packages.speclab.orthopoly1d;
laguerre = opoly.laguerre;
pss = packages.speclab.common.physical_scaleshift_1d;
sss = packages.speclab.common.standard_scaleshift_1d;

opt = laguerre.defaults(varargin{:});
[scale,shift,r] = deal(opt.scale,opt.shift,opt.r);
r = sss(r,opt);

[a,b] = laguerre.coefficients.recurrence(N,opt);
[x,w] = opoly.gauss_radau_quadrature(a,b,N,r);

x = pss(x,opt);
