function[I] = interval(varargin)
% [I] = INTERVAL({SHIFT=0, SCALE=1})
%
%     Returns a 2-element vector specifying the domain of approximation for
%     Jacobi polynomials given the affine scalings SHIFT and SCALE.

global packages;
opt = packages.speclab.orthopoly1d.jacobi.defaults(varargin{:});

I = [-opt.scale+opt.shift, opt.shift+opt.scale];
