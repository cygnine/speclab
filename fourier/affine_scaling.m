function[out] = affine_scaling(interval,varargin)
% [out] = affine_scaling(interval, {N=0, gamma=0, delta=0, resolution_fraction=1})
%
%     Determines the affine parameters scale and shift so that
%     resolution_fraction*N canonical modes lie inside interval. The output
%     out is a struct with two fields, 'scale' and 'shift', specifying the
%     affine map.

global handles;
fourier = handles.speclab.fourier;
inputs = {'N', 'gamma', 'delta', 'resolution_fraction'};
defaults = {0, 0, 0, 1};
opt = handles.common.input_schema(inputs, defaults, [], varargin{:});

out.shift = mean(interval);
if (opt.N==0) | (opt.resolution_fraction==1)
  out.scale = (max(interval) - out.shift)/pi;
  return
end

nodes = fourier.quad.gauss_quadrature(opt.N,'gamma',opt.gamma,'delta',opt.delta);
L = abs(diff(interval))/2;

out.scale = handles.speclab.common.resolution_scaling(L, nodes,...
  'resolution_fraction', opt.resolution_fraction);
