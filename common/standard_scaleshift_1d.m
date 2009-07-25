function[x] = standard_scaleshift_1d(x,varargin)
% [X] = STANDARD_SCALESHIFT1D(X,{SCALE=1,SHIFT=0});
%  
%     Implements the affine shifting+scaling necessary to take things from the
%     one-dimensional physical domain [shift-scale,shift+scale] to the standard
%     domain [-1,1].

global handles;
opt = handles.common.InputSchema({'scale','shift'}, {1,0}, [],varargin{:});

x = (x-opt.shift)/opt.scale;
