function[x,w] = gauss_quadrature(N,varargin)
% gauss_quadrature -- Gauss quadrature rule for Laguerre polynomials
%
% [x,w] = gauss_quadrature(N,{alpha=0,shift=0,scale=1})
%
%     Returns the N-point Gaussian quadrature rule for the Laguerre polynomials. 
%     The weight function is given by speclab.orthopoly1d.laguerre.weights.weight.

global handles;
opoly = handles.speclab.orthopoly1d;
laguerre = opoly.laguerre;
pss = handles.speclab.common.physical_scaleshift_1d;

opt = laguerre.defaults(varargin{:});
[alpha,scale,shift] = deal(opt.alpha,opt.scale,opt.shift);

[a,b] = laguerre.coefficients.recurrence(N+1,opt);
[x,w] = opoly.gauss_quadrature(a,b,N);

x = pss(x,opt);
