function[w] = phase_shifted_sqrt_weight(theta,varargin)
% FUNCTION[W] = PHASE_SHIFTED_SQRT_WEIGHT(THETA,{GAMMA=0, DELTA=0, SCALE=1, SHIFT=0})
%
%     Evaluates the phase-shifted square root weight function for the
%     generalized Fourier Series. The parameters GAMMA and DELTA specify which
%     function class to consider. The affine SHIFT and SCALE parameters indicate
%     the interval over which to evaluate the weight function. Note that this
%     weight function depends on SCALE. 

global handles;
sss = handles.speclab.common.standard_scaleshift_1d;
opt = handles.speclab.fourier.defaults(varargin{:});

theta = sss(theta,opt);

wfourier = handles.speclab.fourier.weights.weight(theta,opt);
w = wfourier.*exp(i*(opt.gamma+opt.delta)/2*(pi-theta));

w = w/sqrt(opt.scale);
