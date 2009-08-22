function[d] = drc_dx(x,varargin)
% [d] = drc_dx(x,varargin)
%
%     The derivative of the mapping x_to_rc. The standard interval rc=[-1,1] is
%     assumed, and the shift and scale parameters are the scaling of the x
%     interval.

global handles;
sss = handles.speclab.common.standard_scaleshift_1d;
opt = handles.common.InputSchema({'shift','scale'}, {0,1}, [],varargin{:});

x = sss(x,opt);

d = 2*(1-x.^2)./(1+x.^2).^2;
d = d/opt.scale;
