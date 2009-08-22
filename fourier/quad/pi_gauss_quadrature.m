function[theta,w] = pi_gauss_quadrature(N,varargin)
% [theta,w] = pi_gauss_quadrature(N,{gamma=0, delta=0, shift=0, scale=1})
%
%     Returns the N-point "pi"-Gauss-Quadrature nodes + weights for the generalized
%     Fourier Series. The default interval is [-pi,pi].
%
%     The parameters gamma, delta determine the weight function for which this
%     quadrature rule is valid. gamma=delta=0 is the canonical Fourier set.
%     shift and scale are affine scaling parameters. 

global handles;
opt = handles.speclab.fourier.defaults(varargin{:});

[theta,w] = fourier.quad.gauss_quadrature(N,opt);

w = w/fourier.weights.weight(theta,opt);
