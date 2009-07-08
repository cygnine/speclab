function[x] = physical_scaleshift1d(x,scale,shift);
% [X] = PHYSICAL_SCALESHIFT1D(X,SCALE,SHIFT);
%  
%     Implements the affine shifting+scaling necessary to take things from the
%     one-dimensional standard domain [-1,1] to the physical domain specified by
%     [shift-scale,shift+scale].

x = x*scale+shift;
