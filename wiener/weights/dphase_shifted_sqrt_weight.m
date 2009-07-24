function[w] = dphase_shifted_sqrt_weight(theta,varargin)
% [W] = DPHASE_SHIFTED_SQRT_WEIGHT(THETA,{S=1, T=0, SCALE=1, SHIFT=0})
%
%     Evaluates the derivative of the phase-shifted square root weight function
%     for the generalized Wiener rational functions. The parameters S and T
%     determine the particular class of functions, and the affine shift factors
%     SCALE and SHIFT determine the affine map of the functions.

global handles;
opt = handles.speclab.wiener.defaults(varargin{:});
psw = handles.speclab.fourier.weights.dphase_shifted_sqrt_weight;
rx = handles.speclab.wiener.maps;

theta = rx.x_to_theta(x,opt);
% Note: it just so happens that dtheta/dx = (1-cos(theta)). Therefore, we
% implicitly take default gamma=s-1 to gamma = (s-1)+1 = s.
fopt = struct('gamma',s,'delta',t);

w = dpsw(theta,fopt);
w = w.*rx.dtheta_dx(x,opt);

% Derivative Jacobian is taken care of in dtheta_dx.
% Weight normalization
% >>>>>> is this right???
w = w/sqrt(opt.scale);
