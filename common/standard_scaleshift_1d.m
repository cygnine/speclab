function[x] = standard_scaleshift1d(x,scale,shift);
% [X] = STANDARD_SCALESHIFT1D(X,SCALE,SHIFT);
%  
%     Implements the affine shifting+scaling necessary to take things from the
%     one-dimensional physical domain [shift-scale,shift+scale] to the standard
%     domain [-1,1].

x = (x-shift)/scale;
