function[x,w] = gauss_radau_quadrature(N,varargin)
% gauss_radau_quadrature -- Gauss-Radau quadrature rule for Laguerre polynomials
%
% [x,w] = gauss_radau_quadrature(N,{alpha=0,shift=0,scale=1,r=0})
%
%     Returns the N-point Gauss-Radau quadrature rule for the Laguerre polynomials.
%     The weight function is given by speclab.orthopoly1d.laguerre.weights.weight.
%     The Radau point is located at x=r.

persistent defaults pss sss recurrence grq
if isempty(defaults)
  from speclab.orthopoly1d.laguerre import defaults
  from speclab.common import standard_scaleshift_1d as sss
  from speclab.common import physical_scaleshift_1d as pss
  from speclab.orthopoly1d.laguerre.coefficients import recurrence
  from speclab.orthopoly1d import gauss_radau_quadrature as grq
end

opt = defaults(varargin{:});
[scale,shift,r] = deal(opt.scale,opt.shift,opt.r);
r = sss(r,opt);

[a,b] = recurrence(N,opt);
[x,w] = grq(a,b,N,r);

x = pss(x,opt);
