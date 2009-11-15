function[x] = standard_scaleshift_1d(x,varargin)
% [x] = standard_scaleshift1d(x,{scale=1,shift=0});
%  
%     Implements the affine shifting+scaling necessary to take things from the
%     one-dimensional physical domain [shift-scale,shift+scale] to the standard
%     domain [-1,1].

persistent strict_inputs
if isempty(strict_inputs)
  from labtools import strict_inputs
end

opt = strict_inputs({'scale','shift'}, {1,0}, [],varargin{:});

x = (x-opt.shift)/opt.scale;
