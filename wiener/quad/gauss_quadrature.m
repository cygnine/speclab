function[x,w] = gauss_quadrature(N,varargin)
% [x,w] = gauss_quadrature(N,{s=1,t=0,scale=1,shift=0})
%
%     Computes an N-point Gauss quadrature rule for the *unweighted* Wiener
%     rational functions. Use pi_gauss_quadrature for the rule for the Wiener
%     rational functions.

global packages;
opt = packages.speclab.wiener.defaults(varargin{:});
rx = packages.speclab.wiener.maps;
fourier = packages.speclab.fourier;

fopt = struct('gamma',opt.s-1, 'delta', opt.t);
% Fourier quadrature rule
[theta,w] = fourier.quad.gauss_quadrature(N,fopt);
x = rx.theta_to_x(theta,opt);
