function[I] = interval(varargin)
% [I] = interval({shift=0, scale=1})
%
%     Returns a 2-element vector specifying the domain of approximation for
%     Laguerre polynomials given the affine scalings shift and scale. 

persistent defaults
if isempty(defaults)
  from speclab.orthopoly1d.laguerre import defaults
end

opt = defaults(varargin{:});

if opt.scale<0
  I = [-Inf, opt.shift];
else
  I = [opt.shift, Inf];
end
