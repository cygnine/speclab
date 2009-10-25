function[w] = dphase_shifted_sqrt_weight(x,varargin)
% [w] = dphase_shifted_sqrt_weight(x,{s=1, t=0, scale=1, shift=0})
%
%     Evaluates the derivative of the phase-shifted square root weight function
%     for the generalized Wiener rational functions. The parameters s and t
%     determine the particular class of functions, and the affine shift factors
%     scale and shift determine the affine map of the functions.

global packages;
opt = packages.speclab.wiener.defaults(varargin{:});
dpsw = packages.speclab.fourier.weights.dphase_shifted_sqrt_weight;
rx = packages.speclab.wiener.maps;

theta = rx.x_to_theta(x,opt);
% Note: it just so happens that dtheta/dx = (1-cos(theta)). Therefore, we
% implicitly take default gamma=s-1 to gamma = (s-1)+1 = s.
fopt = struct('gamma',opt.s,'delta',opt.t);

w = dpsw(theta,fopt);
w = w.*rx.dtheta_dx(x,opt);

% Derivative Jacobian is taken care of in dtheta_dx.
% Weight normalization
% >>>>>> is this right???
w = w/sqrt(opt.scale);
