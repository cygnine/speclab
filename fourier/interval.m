function[I] = interval(varargin)
% [I] = interval({shift=0, scale=1})
%
%     Returns a 2-element vector specifying the domain of approximation for
%     Fourier Series given the affine scalings SHIFT and SCALE. The default
%     interval is [-pi,pi].

persistent defaults
if isempty(defaults)
  from speclab.fourier import defaults
end

opt = defaults(varargin{:});

I = [-pi*opt.scale+opt.shift, opt.shift+pi*opt.scale];
