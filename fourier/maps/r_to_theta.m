function[theta] = r_to_theta(r,varargin)
% function[theta] = r_to_theta(r, {shift=)
%
%     The canonical cos(theta) = r mapping. The standard interval r=[-1,1] is
%     assumed, and the shift and scale parameters refer to the scaling of the
%     theta interval.

global handles;
pss = handles.speclab.common.physical_scaleshift_1d;
opt = handles.common.input_schema({'shift','scale'}, {0,1}, [],varargin{:});

theta = acos(r);
theta = pss(theta,opt);
