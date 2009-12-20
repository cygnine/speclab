function[rc] = x_to_rc(x,varargin)
% [rc] = x_to_rc(x,varargin)
%
%     The canonical (2*x)/(1+x^2) = sin(theta) = rc mapping. The standard
%     interval rc=[-1,1] is assumed, and the shift and scale parameters refer to
%     the scaling of the x interval.

persistent input_schema sss
if isempty(input_schema)
  from labtools import input_schema
  from speclab.common import standard_scaleshift_1d as sss
end

opt = input_schema({'shift','scale'}, {0,1}, [],varargin{:});

x = sss(x,opt);
rc = 2*x./(1+x.^2);
