function[w] = dweight(theta,varargin)
% function[w] = dweight(theta,{gamma=0, delta=0, scale=1, shift=0})
%
%     Evaluates the derivative of the phase-shifted square root weight function
%     for the generalized Fourier Series. The parameters gamma and delta specify
%     which function class to consider. The affine shift and scale parameters
%     indicate the interval over which to evaluate the weight function. Note
%     that this weight function depends on scale. 

persistent dweight dr_dtheta theta_to_r defaults
if isempty(dweight)
  from speclab.fourier import defaults
  from speclab.orthopoly.jacobi.weights import dweight
  from speclab.fourier.maps import dr_dtheta theta_to_r
end

opt = defaults(varargin{:});

dtheta = dr_dtheta(theta,opt);
r = theta_to_r(theta,opt);

w = dweight(r,'alpha',opt.delta, 'beta', opt.gamma);
w = dtheta.*w;
w = w/opt.scale;
