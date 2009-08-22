function[theta] = x_to_theta(x,varargin)
% [theta] = x_to_theta(x,{shift=0,scale=1})
%
%     The canonical x = tan(theta/2) mapping. The standard interval
%     theta=[-pi,pi] is assumed, and the shift and scale parameters refer to the
%     scaling of the x interval.

global handles;
sss = handles.speclab.common.standard_scaleshift_1d;
opt = handles.common.InputSchema({'shift','scale'}, {0,1}, [],varargin{:});

x = sss(x,opt);
theta = 2*atan(x);
