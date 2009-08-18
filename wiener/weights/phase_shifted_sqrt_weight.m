function[w] = phase_shifted_sqrt_weight(x,varargin)
% [w] = phase_shifted_sqrt_weight(x,{s=1, t=0, scale=1, shift=0})
%
%     Evaluates the phase-shifted square root weight function for the
%     generalized Wiener rational functions. The parameters s and t determine
%     the particular class of functions, and the affine shift factors scale and
%     shift determine the affine map of the functions.
%
%     For t=0, scale=1, shift=0, this weight is equivalent to 
%     [i/sqrt(2)*(1 + exp(-i*theta))].^s

global handles;
opt = handles.speclab.wiener.defaults(varargin{:});
psw = handles.speclab.fourier.weights.phase_shifted_sqrt_weight;
rx = handles.speclab.wiener.maps;

theta = rx.x_to_theta(x,opt);
% Note: it just so happens that dtheta/dx = (1-cos(theta)). Therefore, we
% implicitly take default gamma=s-1 to gamma = (s-1)+1 = s.
fopt = struct('gamma',opt.s,'delta',opt.t);

w = psw(theta,fopt);
w = w/sqrt(opt.scale);
