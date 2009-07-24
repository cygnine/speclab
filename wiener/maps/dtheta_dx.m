function[d] = dtheta_dx(x,varargin)
% [D] = DTHETA_DX(X,{SHIFT=0,SCALE=1})
%
%     The derivative of the mapping X_TO_THETA. The standard interval THETA =
%     [-pi,pi] is assumed and SCALE and SHIFT refer to the affine scalings of X. 

global handles;
sss = handles.speclab.common.standard_scaleshift_1d;
opt = handles.common.InputSchema({'shift','scale'}, {0,1}, [],varargin{:});

x = sss(x,opt);
d = 2./(1+x.^2);
d = d/opt.scale;
