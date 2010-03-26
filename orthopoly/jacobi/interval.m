function[I] = interval(varargin)
% [i] = interval({shift=0, scale=1})
%
%     Returns a 2-element vector specifying the domain of approximation for
%     Jacobi polynomials given the affine scalings shift and scale.

persistent defaults
if isempty(defaults)
  from speclab.orthopoly.jacobi import defaults
end

opt = defaults(varargin{:});

I = [-opt.scale+opt.shift, opt.shift+opt.scale];
