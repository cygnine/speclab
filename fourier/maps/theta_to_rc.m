function[r] = theta_to_rc(theta,varargin)
% function[r] = theta_to_r(theta, {shift=0,scale=1})
%
%     The sin(theta) = r mapping (`conjuguate' to r = cos(theta)). The standard
%     interval r=[-1,1] is assumed, and the shift and scale parameters refer to
%     the scaling of the theta interval.
%
%     Note that, in contrast to the mapping r = cos(theta), this mapping depends
%     on the sign parity of theta.

persistent sss input_schema
if isempty(input_schema)
  from labtools import input_schema
  from speclab.common import standard_scaleshift_1d as sss
end

opt = input_schema({'shift','scale'}, {0,1}, [],varargin{:});

theta = sss(theta,opt);
r = sin(theta);
