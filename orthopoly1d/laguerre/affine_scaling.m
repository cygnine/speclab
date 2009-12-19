function[out] = affine_scaling(interval,varargin)
% [out] = affine_scaling(interval, {N=0, alpha=-1/2, resolution_fraction=1})
%
%     Determines the affine parameters scale and shift so that
%     resolution_fraction*N canonical modes lie inside interval. The output
%     out is a struct with two fields, 'scale' and 'shift', specifying the
%     affine map.
%
%     Does not yet support "negative" scaling.

persistent input_schema lag
if isempty(input_schema)
  from labtools import input_schema
  import speclab.orthopoly1d.laguerre as lag
end

%lag = from_as('speclab.orthopoly1d', 'laguerre');

inputs = {'N', 'alpha', 'resolution_fraction'};
defaults = {0, 0, 1};
opt = input_schema(inputs, defaults, [], varargin{:});

out.shift = interval(1);
if (opt.N==0) | (opt.resolution_fraction==1)
  out.scale = 1;
  return
end

%nodes = jac.quad.gauss_quadrature(opt.N,'alpha',opt.alpha,'beta',opt.beta);
%L = abs(diff(interval))/2;

%out.scale = packages.speclab.common.resolution_scaling(L, nodes,...
%  'resolution_fraction', opt.resolution_fraction);
out.scale = 1;
