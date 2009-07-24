function[d] = drc_dx(x,varargin)
% [D] = DRC_DX(X,VARARGIN)
%
%     The derivative of the mapping X_TO_RC. The standard interval rc=[-1,1] is
%     assumed, and the SHIFT and SCALE parameters are the scaling of the X
%     interval.

global handles;
sss = handles.speclab.common.standard_scaleshift_1d;
opt = handles.common.InputSchema({'shift','scale'}, {0,1}, [],varargin{:});

x = sss(x,opt);

d = 2*(1-x.^2)./(1+x.^2).^2;
d = d/opt.scale;
