function[r] = theta_to_r(theta,varargin)
% function[r] = theta_to_r(theta, {shift=0,scale=1})
%
%     The canonical cos(theta) = r mapping. The standard interval r=[-1,1] is
%     assumed, and the shift and scale parameters refer to the scaling of the
%     theta interval.

global handles;
sss = handles.speclab.common.standard_scaleshift_1d;
opt = handles.common.input_schema({'shift','scale'}, {0,1}, [],varargin{:});

theta = sss(theta,opt);
r = cos(theta);
