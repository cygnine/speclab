function[x,w] = gauss_quadrature(N,varargin)
% [x,w] = gauss_quadrature(N,{s=1,t=0,scale=1,shift=0})
%
%     Computes an N-point Gauss quadrature rule for the *unweighted* Wiener
%     rational functions. Use pi_gauss_quadrature for the rule for the Wiener
%     rational functions.

persistent defaults gq theta_to_x
if isempty(defaults)
  from speclab.wiener import defaults
  from speclab.wiener.maps import theta_to_x
  from speclab.fourier.quad import gauss_quadrature as gq
end

opt = defaults(varargin{:});

fopt = struct('gamma',opt.s-1, 'delta', opt.t);
% Fourier quadrature rule
[theta,w] = gq(N,fopt);
x = theta_to_x(theta,opt);
