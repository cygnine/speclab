function[x,w] = gauss_quadrature(N,varargin)
% gauss_quadrature -- Gauss quadrature rule for Hermite polynomials
%
% [x,w] = gauss_quadrature(N,{mu=0,shift=0,scale=1})
%
%     Returns the N-point Gaussian quadrature rule for the Hermite polynomials. 
%     The weight function is given by speclab.orthopoly1d.hermite.weights.weight.

global handles;
opoly = handles.speclab.orthopoly1d;
hermite = opoly.hermite;
pss = handles.speclab.common.physical_scaleshift_1d;

opt = hermite.defaults(varargin{:});
[mu,scale,shift] = deal(opt.mu,opt.scale,opt.shift);

[a,b] = hermite.coefficients.recurrence(N+1,opt);
[x,w] = opoly.gauss_quadrature(a,b,N);

x = pss(x,opt);
