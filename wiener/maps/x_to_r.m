function[r] = x_to_r(x,varargin)
% [R] = X_TO_R(X,VARARGIN)
%
%     The canonical (1-x^2)/(1+x^2) = cos(theta) = r mapping. The standard
%     interval r=[-1,1] is assumed, and the SHIFT and SCALE parameters refer to
%     the scaling of the X interval.

global handles;
sss = handles.speclab.common.standard_scaleshift_1d;
opt = handles.common.InputSchema({'shift','scale'}, {0,1}, [],varargin{:});

x = sss(x,opt);
r = (1-x.^2)./(1+x.^2);
