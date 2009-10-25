function[d] = drc_dtheta(theta,varargin)
% function[d] = drc_dtheta(theta, {shift=0,scale=1})
%
%     The Jacobian of the sin(theta) = r mapping (`conjuguate' to r =
%     cos(theta)). The standard interval r=[-1,1] is assumed, and the shift and
%     scale parameters refer to the scaling of the theta interval.
%
%     Note that, in contrast to the mapping r = cos(theta), this mapping depends
%     on the sign parity of theta.

global packages;
sss = packages.speclab.common.standard_scaleshift_1d;
opt = packages.labtools.input_schema({'shift','scale'}, {0,1}, [],varargin{:});

theta = sss(theta,opt);
d = cos(theta)/opt.scale;
