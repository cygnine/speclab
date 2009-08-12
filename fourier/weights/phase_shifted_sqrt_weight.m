function[w] = phase_shifted_sqrt_weight(theta,varargin)
% function[w] = phase_shifted_sqrt_weight(theta,{gamma=0, delta=0, scale=1, shift=0})
%
%     Evaluates the phase-shifted square root weight function for the
%     generalized Fourier Series. The parameters gamma and delta specify which
%     function class to consider. The affine shift and scale parameters indicate
%     the interval over which to evaluate the weight function. Note that this
%     weight function depends on scale. 

global handles;
sss = handles.speclab.common.standard_scaleshift_1d;
opt = handles.speclab.fourier.defaults(varargin{:});

opt.gamma = opt.gamma/2;
opt.delta = opt.delta/2;
wfourier = handles.speclab.fourier.weights.weight(theta,opt);
theta = sss(theta,opt);
w = wfourier.*exp(i*(opt.gamma+opt.delta)*(pi-theta));

w = w/sqrt(opt.scale);
