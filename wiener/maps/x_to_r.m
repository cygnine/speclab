function[r] = x_to_r(x,varargin)
% [r] = x_to_r(x,varargin)
%
%     The canonical (1-x^2)/(1+x^2) = cos(theta) = r mapping. The standard
%     interval r=[-1,1] is assumed, and the shift and scale parameters refer to
%     the scaling of the X interval.

global packages;
sss = packages.speclab.common.standard_scaleshift_1d;
opt = packages.common.input_schema({'shift','scale'}, {0,1}, [],varargin{:});

x = sss(x,opt);
r = (1-x.^2)./(1+x.^2);
