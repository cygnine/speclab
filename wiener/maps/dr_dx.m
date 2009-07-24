function[d] = dr_dx(x,varargin)
% [D] = DR_DX(X,VARARGIN)
%
%     The derivative of the mapping X_TO_R. The standard interval r=[-1,1] is
%     assumed, and the SHIFT and SCALE parameters refer to the scaling of the X
%     interval.

global handles;
sss = handles.speclab.common.standard_scaleshift_1d;
opt = handles.common.InputSchema({'shift','scale'}, {0,1}, [],varargin{:});

x = sss(x,opt);
d = -4*x./(x.^2+1).^2;
d = d/opt.scale;
