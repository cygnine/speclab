function[out] = affine_scaling(interval,varargin)
% [out] = affine_scaling(interval, {N=0, alpha=-1/2, beta=-1/2, resolution_fraction=1})
%
%     Determines the affine parameters scale and shift so that
%     resolution_fraction*N canonical modes lie inside interval. The output
%     out is a struct with two fields, 'scale' and 'shift', specifying the
%     affine map.

persistent input_schema gauss_qudarature resolution_scaling
if isempty(input_schema)
  from labtools import input_schema
  from speclab.orthopoly.jacobi.quad import gauss_quadrature
  from speclab.common import resolution_scaling
end

inputs = {'N', 'alpha', 'beta', 'resolution_fraction'};
defaults = {0, -1/2, -1/2, 1};
opt = input_schema(inputs, defaults, [], varargin{:});

out.shift = mean(interval);
if (opt.N==0) | (opt.resolution_fraction==1)
  out.scale = max(interval) - out.shift;
  return
end

nodes = gauss_quadrature(opt.N,'alpha',opt.alpha,'beta',opt.beta);
L = abs(diff(interval))/2;

out.scale = resolution_scaling(L, nodes, 'resolution_fraction', opt.resolution_fraction);
