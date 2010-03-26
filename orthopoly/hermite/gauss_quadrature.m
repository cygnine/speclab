function[x,w] = gauss_quadrature(N,varargin)
% gauss_quadrature -- Gauss quadrature rule for Hermite polynomials
%
% [x,w] = gauss_quadrature(N,{mu=0,shift=0,scale=1})
%
%     Returns the N-point Gaussian quadrature rule for the Hermite polynomials. 
%     The weight function is given by speclab.orthopoly.hermite.weights.weight.

persistent recurrence gq pss defaults
if isempty(recurrence)
  from speclab.orthopoly.hermite import defaults recurrence
  from speclab.orthopoly import gauss_quadrature as gq
  from speclab.common import physical_scaleshift_1d as pss
end

opt = defaults(varargin{:});

[a,b] = recurrence(N+1,opt);
[x,w] = gq(a,b,N);

switch opt.weight_normalization
case 'probability'
  w = w/b(1);
end

x = pss(x,opt);
