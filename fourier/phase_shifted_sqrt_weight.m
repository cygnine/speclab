function[w] = phase_shifted_sqrt_weight(theta,varargin)
% function[w] = phase_shifted_sqrt_weight(theta,{gamma=0, delta=0, scale=1, shift=0})
%
%     Evaluates the phase-shifted square root weight function for the
%     generalized Fourier Series. The parameters gamma and delta specify which
%     function class to consider. The affine shift and scale parameters indicate
%     the interval over which to evaluate the weight function. Note that this
%     weight function depends on scale. 

persistent sss defaults weight
if isempty(defaults)
  from speclab.common import standard_scaleshift_1d as sss
  from speclab.fourier import defaults
  from speclab.fourier import weight
end

opt = defaults(varargin{:});

opt.gamma = opt.gamma/2;
opt.delta = opt.delta/2;
wfourier = weight(theta,opt);
theta = sss(theta,opt);
w = wfourier.*exp(i*(opt.gamma+opt.delta)*(pi-theta));

w = w/sqrt(opt.scale);
