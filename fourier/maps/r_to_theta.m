function[theta] = r_to_theta(r,varargin)
% FUNCTION[THETA] = R_TO_THETA(R, {shift=)
%
%     The canonical cos(theta) = r mapping. The standard interval r=[-1,1] is
%     assumed, and the SHIFT and SCALE parameters refer to the scaling of the
%     THETA interval.

global handles;
pss = handles.speclab.common.physical_scaleshift_1d;
opt = handles.common.InputSchema({'shift','scale'}, {0,1}, [],varargin{:});

theta = acos(r);
theta = pss(theta,opt);
