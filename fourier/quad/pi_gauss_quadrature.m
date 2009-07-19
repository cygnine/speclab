function[theta,w] = pi_gauss_quadrature(N,varargin)
% [THETA,W] = PI_GAUSS_QUADRATURE(N,{GAMMA=0, DELTA=0, SHIFT=0, SCALE=1})
%
%     Returns the "pi"-Gauss-Quadrature nodes + weights for the generalized
%     Fourier Series. The default interval is [0,2*pi].
%
%     The parameters GAMMA, DELTA determine the weight function for which this
%     quadrature rule is valid. GAMMA=DELTA=0 is the canonical Fourier set.
%     SHIFT and SCALE are affine scaling parameters. 

global handles;
opt = handles.speclab.fourier.defaults(varargin{:});

[theta,w] = fourier.quad.gauss_quadrature(N,opt);

w = w/fourier.weights.weight(theta,opt);
