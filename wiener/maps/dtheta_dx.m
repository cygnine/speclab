function[d] = dtheta_dx(x,varargin)
% [d] = dtheta_dx(x,{shift=0,scale=1})
%
%     The derivative of the mapping x_to_theta. The standard interval theta =
%     [-pi,pi] is assumed and scale and shift refer to the affine scalings of x. 

persistent input_schema sss
if isempty(input_schema)
  from labtools import input_schema
  from speclab.common import standard_scaleshift_1d as sss
end

opt = input_schema({'shift','scale'}, {0,1}, [],varargin{:});

x = sss(x,opt);
d = 2./(1+x.^2);
d = d/opt.scale;
