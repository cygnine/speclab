function[x] = physical_scaleshift(x,varargin)
% [x] = physical_scaleshift(x,{scale=1,shift=0});
%  
%     Implements the affine shifting+scaling necessary to take things from the
%     rectangualr physical domain [shift-scale,shift+scale] to the standard
%     domain [-1,1]. Supports multi-dimensionality.
%
%     Superscedes physical_scaleshift_1d, which will soon be deprecated.

persistent strict_inputs spdiag
if isempty(strict_inputs)
  from labtools import strict_inputs spdiag
end

opt = strict_inputs({'scale','shift'}, {1,0}, [],varargin{:});

dim = length(opt.scale);
assert(length(opt.shift)==dim, 'The shift and scale must have the same lengths');
assert(size(x,2)==dim, 'The input points must have the correct dimensionality');

if dim==1
  xsize = size(x);
  x = x(:);
end

N = size(x,1);

x = x*spdiag(opt.scale(:)) + repmat(opt.shift(:).', [N 1]);

if dim==1
  x = reshape(x, xsize);
end
