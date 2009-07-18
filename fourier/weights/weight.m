function[w] = weight(theta,varargin)
% FUNCTION[W] = WEIGHT(THETA,{GAMMA=0, DELTA=0, SCALE=1, SHIFT=0})
%
%     Evaluates the weight function for the generalized Fourier Series. The
%     parameters GAMMA and DELTA specify which function class to consider. The
%     affine SHIFT and SCALE parameters indicate the interval over which to
%     evaluate the weight function. Note that this weight function depends on
%     SCALE. 

global handles;
sss = handles.speclab.common.standard_scaleshift_1d;
opt = handles.speclab.fourier.defaults(varargin{:});

theta = sss(theta,opt);

w = ((1-cos(theta)).^opt.gamma) .* ((1+cos(theta)).^opt.delta);
w = w/opt.scale;
