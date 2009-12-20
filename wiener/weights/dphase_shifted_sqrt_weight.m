function[w] = dphase_shifted_sqrt_weight(x,varargin)
% [w] = dphase_shifted_sqrt_weight(x,{s=1, t=0, scale=1, shift=0})
%
%     Evaluates the derivative of the phase-shifted square root weight function
%     for the generalized Wiener rational functions. The parameters s and t
%     determine the particular class of functions, and the affine shift factors
%     scale and shift determine the affine map of the functions.

persistent defaults dsqrt_weight x_to_theta dtheta_dx
if isempty(defaults)
  from speclab.wiener import defaults
  from speclab.fourier.weights import dphase_shifted_sqrt_weight as dsqrt_weight
  from speclab.wiener.maps import x_to_theta dtheta_dx
end

opt = defaults(varargin{:});

theta = x_to_theta(x,opt);
% Note: it just so happens that dtheta/dx = (1-cos(theta)). Therefore, we
% implicitly take default gamma=s-1 to gamma = (s-1)+1 = s.
fopt = struct('gamma',opt.s,'delta',opt.t);

w = dsqrt_weight(theta,fopt);
w = w.*dtheta_dx(x,opt);

% Derivative Jacobian is taken care of in dtheta_dx.
% Weight normalization
% >>>>>> is this right???
w = w/sqrt(opt.scale);
