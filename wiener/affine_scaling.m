function[out] = affine_scaling(interval,N,varargin)
% [out] = affine_scaling(interval, N {s=1, t=0, resolution_fraction=1})
%
%     Determines the affine parameters scale and shift so that
%     resolution_fraction*N canonical modes lie inside interval. The output
%     out is a struct with two fields, 'scale' and 'shift', specifying the
%     affine map.

global packages;
wiener = packages.speclab.wiener;
inputs = {'s', 't', 'resolution_fraction'};
defaults = {1, 0, 1};
opt = packages.common.input_schema(inputs, defaults, [], varargin{:});

out.shift = mean(interval);
if (N<1)
  error('The number of degrees of freedom must be > 0');
  return
end

nodes = wiener.quad.gauss_quadrature(N,'s',opt.s,'t',opt.t);
L = abs(diff(interval))/2;

out.scale = packages.speclab.common.resolution_scaling(L, nodes,...
  'resolution_fraction', opt.resolution_fraction);
