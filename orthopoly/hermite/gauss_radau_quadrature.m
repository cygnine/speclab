function[x,w] = gauss_radau_quadrature(N,varargin)
% gauss_radau_quadrature -- Gauss-Radau quadrature rule for Hermite polynomials
%
% [x,w] = gauss_radau_quadrature(N,{mu=0,shift=0,scale=1,r=0})
%
%     Returns the N-point Gauss-Radau quadrature rule for the Hermite polynomials.
%     The weight function is given by speclab.orthopoly.hermite.weights.weight.
%     The Radau point is located at x=r.

persistent defaults pss sss recurrence grq
if isempty(defaults)
  from speclab.orthopoly.hermite import defaults recurrence
  from speclab.common import physical_scaleshift_1d as pss
  from speclab.common import standard_scaleshift_1d as sss
  from speclab.orthopoly import gauss_radau_quadrature as grq
end

opt = defaults(varargin{:});
opt.x = sss(opt.x,opt);

[a,b] = recurrence(N,opt);
[x,w] = grq(a,b,N,opt.x);

switch opt.weight_normalization
case 'probability'
  w = w/b(1);
end

x = pss(x,opt);
