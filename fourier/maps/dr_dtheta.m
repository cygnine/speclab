function[d] = dr_dtheta(theta,varargin)
% function[d] = dr_dtheta(theta, {shift=0,scale=1})
%
%     The Jacobian of the canonical cos(theta) = r mapping. The standard
%     interval r=[-1,1] is assumed, and the shift and scale parameters refer to
%     the scaling of the theta interval.

persistent sss input_schema
if isempty(input_schema)
  from labtools import input_schema
  from speclab.common import standard_scaleshift_1d as sss
end

opt = input_schema({'shift','scale'}, {0,1}, [],varargin{:});

theta = sss(theta,opt);
d = -sin(theta)/opt.scale;
