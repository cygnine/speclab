function[r] = theta_to_rc(theta,varargin)
% FUNCTION[R] = THETA_TO_R(THETA, {SHIFT=0,SCALE=1})
%
%     The sin(theta) = r mapping (`conjuguate' to r = cos(theta)). The standard
%     interval r=[-1,1] is assumed, and the SHIFT and SCALE parameters refer to
%     the scaling of the THETA interval.
%
%     Note that, in contrast to the mapping r = cos(theta), this mapping depends
%     on the sign parity of THETA.

global handles;
sss = handles.speclab.common.standard_scaleshift_1d;
opt = handles.common.InputSchema({'shift','scale'}, {0,1}, [],varargin{:});

theta = sss(theta,opt);
r = sin(theta);
