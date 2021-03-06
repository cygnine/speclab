function[d] = dr_dx(x,varargin)
% [d] = dr_dx(x,varargin)
%
%     The derivative of the mapping x_to_r. The standard interval r=[-1,1] is
%     assumed, and the shift and scale parameters refer to the scaling of the x
%     interval.

persistent input_schema sss
if isempty(input_schema)
  from labtools import input_schema
  from speclab.common import standard_scaleshift_1d as sss
end

opt = input_schema({'shift','scale'}, {0,1}, [],varargin{:});

x = sss(x,opt);
d = -4*x./(x.^2+1).^2;
d = d/opt.scale;
