function[theta] = x_to_theta(x,varargin)
% [theta] = x_to_theta(x,{shift=0,scale=1})
%
%     The canonical x = tan(theta/2) mapping. The standard interval
%     theta=[-pi,pi] is assumed, and the shift and scale parameters refer to the
%     scaling of the x interval.

persistent input_schema sss
if isempty(input_schema)
  from labtools import input_schema
  from speclab.common import standard_scaleshift_1d as sss
end

opt = input_schema({'shift','scale'}, {0,1}, [],varargin{:});

x = sss(x,opt);
theta = 2*atan(x);
