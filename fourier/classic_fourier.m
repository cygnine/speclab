function[tf] = classic_fourier(varargin)
% [tf] = classic_fourier({gamma=0, delta=0})
%
%     Determines whether the parameters gamma and delta are 0, which means that
%     the corresponding Fourier Series family is the canonical one.

global handles;
opt = handles.common.InputSchema({'gamma', 'delta'}, {0,0}, [], varargin{:});

tol = 1e-12;
tf = abs(opt.gamma)<tol & abs(opt.delta)<tol;
