function[x] = physical_scaleshift_1d(x,varargin)
% [x] = physical_scaleshift1d(x,{scale=1,shift=0});
%  
%     Implements the affine shifting+scaling necessary to take things from the
%     one-dimensional standard domain [-1,1] to the physical domain specified by
%     [shift-scale,shift+scale].

persistent strict_inputs
if isempty(strict_inputs)
  from labtools import strict_inputs;
end

opt = strict_inputs({'scale','shift'}, {1,0}, [],varargin{:});

x = x*opt.scale+opt.shift;
