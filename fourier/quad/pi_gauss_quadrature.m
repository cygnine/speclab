function[theta,w] = pi_gauss_quadrature(N,varargin)
% [theta,w] = pi_gauss_quadrature(N,{gamma=0, delta=0, shift=0, scale=1})
%
%     Returns the N-point "pi"-Gauss-Quadrature nodes + weights for the generalized
%     Fourier Series. The default interval is [-pi,pi].
%
%     The parameters gamma, delta determine the weight function for which this
%     quadrature rule is valid. gamma=delta=0 is the canonical Fourier set.
%     shift and scale are affine scaling parameters. 

persistent gq defaults weight
if isempty(gq)
  from speclab.fourier import defaults
  from speclab.fourier.quad import gauss_quadrature as gq
  from speclab.fourier.weights import weight
end

opt = defaults(varargin{:});

[theta,w] = gq(N,opt);

w = w/weight(theta,opt);
