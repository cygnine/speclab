function[d] = dtheta_dr(r,varargin)
% FUNCTION[D] = DTHETA_DR(R, {SHIFT=0,SCALE=1})
%
%     The Jacobian of the canonical cos(theta) = r mapping. The standard
%     interval r=[-1,1] is assumed, and the SHIFT and SCALE parameters refer to
%     the scaling of the THETA interval.

global handles;
opt = handles.common.InputSchema({'shift','scale'}, {0,1}, [],varargin{:});

d = -opt.scale/(sqrt(1-r.^2));
