function[r] = theta_to_r(theta,varargin)
% FUNCTION[R] = THETA_TO_R(THETA, {SHIFT=0,SCALE=1})
%
%     The canonical cos(theta) = r mapping. The standard interval r=[-1,1] is
%     assumed, and the SHIFT and SCALE parameters refer to the scaling of the
%     THETA interval.

global handles;
sss = handles.speclab.common.standard_scaleshift_1d;
opt = handles.common.InputSchema({'shift','scale'}, {0,1}, [],varargin{:});

theta = sss(theta,opt);
r = cos(theta);
