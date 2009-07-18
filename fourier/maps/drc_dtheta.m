function[d] = drc_dtheta(theta,varargin)
% FUNCTION[D] = DRC_DTHETA(THETA, {SHIFT=0,SCALE=1})
%
%     The Jacobian of the sin(theta) = r mapping (`conjuguate' to r =
%     cos(theta)). The standard interval r=[-1,1] is assumed, and the SHIFT and
%     SCALE parameters refer to the scaling of the THETA interval.
%
%     Note that, in contrast to the mapping r = cos(theta), this mapping depends
%     on the sign parity of THETA.

global handles;
sss = handles.speclab.common.standard_scaleshift_1d;
opt = handles.common.InputSchema({'shift','scale'}, {0,1}, [],varargin{:});

theta = sss(theta,opt);
d = cos(theta)*opt.scale;
