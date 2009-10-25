function[d] = dtheta_dr(r,varargin)
% function[d] = dtheta_dr(r, {shift=0,scale=1})
%
%     The Jacobian of the canonical cos(theta) = r mapping. The standard
%     interval r=[-1,1] is assumed, and the shift and scale parameters refer to
%     the scaling of the theta interval.

global packages;
opt = packages.labtools.input_schema({'shift','scale'}, {0,1}, [],varargin{:});

d = -opt.scale/(sqrt(1-r.^2));
