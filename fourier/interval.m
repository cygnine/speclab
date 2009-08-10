function[I] = interval(varargin)
% [I] = INTERVAL({SHIFT=0, SCALE=1})
%
%     Returns a 2-element vector specifying the domain of approximation for
%     Fourier Series given the affine scalings SHIFT and SCALE. The default
%     interval is [-pi,pi].

global handles;
opt = handles.speclab.fourier.defaults(varargin{:});

I = [-pi*opt.scale+opt.shift, opt.shift+pi*opt.scale];
