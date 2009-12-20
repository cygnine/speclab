function[x,w] = gauss_quadrature(N,varargin)
% gauss_quadrature -- Gauss quadrature rule for Laguerre polynomials
%
% [x,w] = gauss_quadrature(N,{alpha=0,shift=0,scale=1})
%
%     Returns the N-point Gaussian quadrature rule for the Laguerre polynomials. 
%     The weight function is given by speclab.orthopoly1d.laguerre.weights.weight.

persistent defaults recurrence gq pss
if isempty(defaults)
  from speclab.orthopoly1d.laguerre import defaults
  from speclab.orthopoly1d.laguerre.coefficients import recurrence
  from speclab.orthopoly1d import gauss_quadrature as gq
  from speclab.common import physical_scaleshift_1d as pss
end

opt = defaults(varargin{:});
[alpha,scale,shift] = deal(opt.alpha,opt.scale,opt.shift);

[a,b] = recurrence(N+1,opt);
[x,w] = gq(a,b,N);

x = pss(x,opt);
