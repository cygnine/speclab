function[d] = dr_dtheta(theta,varargin)
% FUNCTION[D] = DR_DTHETA(THETA, {SHIFT=0,SCALE=1})
%
%     The Jacobian of the canonical cos(theta) = r mapping. The standard
%     interval r=[-1,1] is assumed, and the SHIFT and SCALE parameters refer to
%     the scaling of the THETA interval.

global handles;
sss = handles.speclab.common.standard_scaleshift_1d;
opt = handles.common.InputSchema({'shift','scale'}, {0,1}, [],varargin{:});

theta = sss(theta,opt);
d = -sin(theta)/opt.scale;
