function[d] = dr_dx(x,varargin)
% [d] = dr_dx(x,varargin)
%
%     The derivative of the mapping x_to_r. The standard interval r=[-1,1] is
%     assumed, and the shift and scale parameters refer to the scaling of the x
%     interval.

global handles;
sss = handles.speclab.common.standard_scaleshift_1d;
opt = handles.common.InputSchema({'shift','scale'}, {0,1}, [],varargin{:});

x = sss(x,opt);
d = -4*x./(x.^2+1).^2;
d = d/opt.scale;
