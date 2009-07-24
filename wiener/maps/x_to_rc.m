function[rc] = x_to_rc(x,varargin)
% [RC] = X_TO_RC(X,VARARGIN)
%
%     The canonical (2*x)/(1+x^2) = sin(theta) = rc mapping. The standard
%     interval rc=[-1,1] is assumed, and the SHIFT and SCALE parameters refer to
%     the scaling of the X interval.

global handles;
sss = handles.speclab.common.standard_scaleshift_1d;
opt = handles.common.InputSchema({'shift','scale'}, {0,1}, [],varargin{:});

x = sss(x,opt);
rc = 2*x./(1+x.^2);
