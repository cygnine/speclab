function[x,w] = gauss_radau_quadrature(N,varargin)
% gauss_radau_quadrature -- Gauss-Radau quadrature rule for Hermite polynomials
%
% [x,w] = gauss_radau_quadrature(N,{mu=0,shift=0,scale=1,r=0})
%
%     Returns the N-point Gauss-Radau quadrature rule for the Hermite polynomials.
%     The weight function is given by speclab.orthopoly1d.hermite.weights.weight.
%     The Radau point is located at x=r.

global packages;
opoly = packages.speclab.orthopoly1d;
hermite = opoly.hermite;
pss = packages.speclab.common.physical_scaleshift_1d;
sss = packages.speclab.common.standard_scaleshift_1d;

opt = hermite.defaults(varargin{:});
[mu,scale,shift,r] = deal(opt.scale,opt.shift,opt.x);
r = sss(r,opt);

[a,b] = hermite.coefficients.recurrence(N,opt);
[x,w] = opoly.gauss_radau_quadrature(a,b,N,r);

x = pss(x,opt);
