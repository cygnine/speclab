function[x] = theta_to_x(theta,varargin)
% [X] = THETA_TO_X(THETA,{SHIFT=0,SCALE=1})
%
%     The canonical x = tan(theta/2) mapping. The standard interval
%     theta=[-pi,pi] is assumed, and the SHIFT and SCALE parameters refer to the
%     scaling of the X interval.

global handles;
pss = handles.speclab.common.physical_scaleshift_1d;
opt = handles.common.InputSchema({'shift','scale'}, {0,1}, [],varargin{:});

x = tan(theta/2);
x = pss(x,opt);
