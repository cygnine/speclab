function[out] = affine_scaling(interval,varargin)
% [out] = affine_scaling(interval, {N=0, alpha=-1/2, beta=-1/2, resolution_fraction=1})
%
%     Determines the affine parameters scale and shift so that
%     resolution_fraction*N canonical modes lie inside interval. The output
%     out is a struct with two fields, 'scale' and 'shift', specifying the
%     affine map.

global packages;
jac = packages.speclab.orthopoly1d.jacobi;
inputs = {'N', 'alpha', 'beta', 'resolution_fraction'};
defaults = {0, -1/2, -1/2, 1};
opt = packages.labtools.input_schema(inputs, defaults, [], varargin{:});

out.shift = mean(interval);
if (opt.N==0) | (opt.resolution_fraction==1)
  out.scale = max(interval) - out.shift;
  return
end

nodes = jac.quad.gauss_quadrature(opt.N,'alpha',opt.alpha,'beta',opt.beta);
L = abs(diff(interval))/2;

out.scale = packages.speclab.common.resolution_scaling(L, nodes,...
  'resolution_fraction', opt.resolution_fraction);
