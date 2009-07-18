function[w] = dphase_shifted_sqrt_weight(theta,varargin)
% FUNCTION[W] = DPHASE_SHIFTED_SQRT_WEIGHT(THETA,{GAMMA=0, DELTA=0, SCALE=1, SHIFT=0})
%
%     Evaluates the derivative of the phase-shifted square root weight function
%     for the generalized Fourier Series. The parameters GAMMA and DELTA specify
%     which function class to consider. The affine SHIFT and SCALE parameters
%     indicate the interval over which to evaluate the weight function. Note
%     that this weight function depends on SCALE. 

global handles;
sss = handles.speclab.common.standard_scaleshift_1d;
opt = handles.speclab.fourier.defaults(varargin{:});

theta = sss(theta,opt);
tol = 1e-6;

phase = exp(i*(opt.gamma+opt.delta)/2*(pi-theta))*...
        2^((opt.gamma+opt.delta-4)/2);

if abs(opt.delta)>tol
  phase = phase * abs(sin(theta/2).^(opt.delta-1));
end
if abs(opt.gamma)>tol
  phase = phase * cos(theta/2).^(opt.gamma-1);
end

w = phase * (opt.delta*(1+exp(-i*theta)) - opt.gamma*(1+exp(i*theta)));

% Weight normalization + derivative Jacobian
w = w/opt.scale.^2;
