function[w] = phase_shifted_sqrt_weight(theta,varargin)
% [W] = PHASE_SHIFTED_SQRT_WEIGHT(THETA,{S=1, T=0, SCALE=1, SHIFT=0})
%
%     Evaluates the phase-shifted square root weight function for the
%     generalized Wiener rational functions. The parameters S and T determine
%     the particular class of functions, and the affine shift factors SCALE and
%     SHIFT determine the affine map of the functions.

global handles;
opt = handles.speclab.wiener.defaults(varargin{:});
psw = handles.speclab.fourier.weights.phase_shifted_sqrt_weight;
rx = handles.speclab.wiener.maps;

theta = rx.x_to_theta(x,opt);
% Note: it just so happens that dtheta/dx = (1-cos(theta)). Therefore, we
% implicitly take default gamma=s-1 to gamma = (s-1)+1 = s.
fopt = struct('gamma',s,'delta',t);

w = psw(theta,fopt);
w = w/sqrt(opt.scale);
