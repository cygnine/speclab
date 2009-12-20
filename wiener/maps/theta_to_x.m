function[x] = theta_to_x(theta,varargin)
% [x] = theta_to_x(theta,{shift=0,scale=1})
%
%     The canonical x = tan(theta/2) mapping. The standard interval
%     theta=[-pi,pi] is assumed, and the shift and scale parameters refer to the
%     scaling of the x interval.

persistent input_schema pss
if isempty(input_schema)
  from labtools import input_schema
  from speclab.common import physical_scaleshift_1d as pss
end

opt = input_schema({'shift','scale'}, {0,1}, [],varargin{:});

x = tan(theta/2);
x = pss(x,opt);
