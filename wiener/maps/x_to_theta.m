function[theta] = x_to_theta(x,varargin)
% [THETA] = X_TO_THETA(X,{SHIFT=0,SCALE=1})
%
%     The canonical x = tan(theta/2) mapping. The standard interval
%     theta=[-pi,pi] is assumed, and the SHIFT and SCALE parameters refer to the
%     scaling of the X interval.

global handles;
sss = handles.speclab.common.standard_scaleshift_1d;
opt = handles.common.InputSchema({'shift','scale'}, {0,1}, [],varargin{:});

x = sss(x,opt);
theta = 2*atan(x);
