function[rc] = x_to_rc(x,varargin)
% [rc] = x_to_rc(x,varargin)
%
%     The canonical (2*x)/(1+x^2) = sin(theta) = rc mapping. The standard
%     interval rc=[-1,1] is assumed, and the shift and scale parameters refer to
%     the scaling of the x interval.

global handles;
sss = handles.speclab.common.standard_scaleshift_1d;
opt = handles.common.input_schema({'shift','scale'}, {0,1}, [],varargin{:});

x = sss(x,opt);
rc = 2*x./(1+x.^2);
