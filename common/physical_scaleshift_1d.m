function[x] = physical_scaleshift_1d(x,varargin)
% [X] = PHYSICAL_SCALESHIFT1D(X,{SCALE=1,SHIFT=0});
%  
%     Implements the affine shifting+scaling necessary to take things from the
%     one-dimensional standard domain [-1,1] to the physical domain specified by
%     [shift-scale,shift+scale].

global handles;
opt = handles.common.input_schema({'scale','shift'}, {1,0}, [],varargin{:});

x = x*opt.scale+opt.shift;
