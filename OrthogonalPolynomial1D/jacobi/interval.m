function[I] = interval(varargin)
% [I] = INTERVAL({SHIFT=0, SCALE=1})
%
%     Returns a 2-element vector specifying the domain of approximation for
%     Jacobi polynomials given the affine scalings SHIFT and SCALE.

global handles;
opt = handles.speclab.OrthogonalPolynomial1D.jacobi.defaults(varargin{:});

I = [-opt.scale+opt.shift, opt.shift+opt.scale];
