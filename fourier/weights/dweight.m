function[w] = dweight(theta,varargin)
% function[w] = dweight(theta,{gamma=0, delta=0, scale=1, shift=0})
%
%     Evaluates the derivative of the phase-shifted square root weight function
%     for the generalized Fourier Series. The parameters gamma and delta specify
%     which function class to consider. The affine shift and scale parameters
%     indicate the interval over which to evaluate the weight function. Note
%     that this weight function depends on scale. 

global packages;
opt = packages.speclab.fourier.defaults(varargin{:});
jac = packages.speclab.orthopoly1d.jacobi;
tr = packages.speclab.fourier.maps;

dtheta = tr.dr_dtheta(theta,opt);
r = tr.theta_to_r(theta,opt);

w = jac.weights.dweight(r,'alpha',opt.delta, 'beta', opt.gamma);
w = dtheta.*w;
w = w/opt.scale;
