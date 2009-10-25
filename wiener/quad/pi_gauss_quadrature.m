function[x,w] = pi_gauss_quadrature(N,varargin)
% [x,w] = pi_gauss_quadrature(N,{s=1,t=0,scale=1,shift=0})
%
%     Computes the N-point pi-Gauss quadrature rule for the Wiener
%     rational functions. 

global packages;
opt = packages.speclab.wiener.defaults(varargin{:});
wiener = packages.speclab.wiener;

% Unweighted quadrature rule
[x,w] = wiener.quad.gauss_quadrature(N,opt);

weight_factor = wiener.weights.weight(x,opt);
w = w./weight_factor;
