function[x,w] = pi_gauss_quadrature(N,varargin)
% [x,w] = pi_gauss_quadrature(N,{s=1,t=0,scale=1,shift=0})
%
%     Computes the N-point pi-Gauss quadrature rule for the Wiener
%     rational functions. 

persistent defaults gauss_quadrature weight
if isempty(defaults)
  from speclab.wiener import defaults
  from speclab.wiener.quad import gauss_quadrature
  from speclab.wiener.weights import weight
end

opt = defaults(varargin{:});

% Unweighted quadrature rule
[x,w] = gauss_quadrature(N,opt);

weight_factor = weight(x,opt);
w = w./weight_factor;
