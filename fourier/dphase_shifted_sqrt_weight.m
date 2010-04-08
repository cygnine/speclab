function[w] = dphase_shifted_sqrt_weight(theta,varargin)
% function[w] = dphase_shifted_sqrt_weight(theta,{gamma=0, delta=0, scale=1, shift=0})
%
%     Evaluates the derivative of the phase-shifted square root weight function
%     for the generalized Fourier Series. The parameters gamma and delta specify
%     which function class to consider. The affine shift and scale parameters
%     indicate the interval over which to evaluate the weight function. Note
%     that this weight function depends on scale. 

persistent sss weight dweight defaults
if isempty(weight)
  from speclab.fourier import defaults
  from speclab.common import standard_scaleshift_1d as sss
  from speclab.fourier import weight dweight
end

opt = defaults(varargin{:});
%tr = packages.speclab.fourier.maps;

opt.gamma = opt.gamma/2;
opt.delta = opt.delta/2;
fw = weight(theta,opt);
dfw = dweight(theta,opt);

theta = sss(theta,opt);

phase = exp(i*(opt.gamma+opt.delta)*(pi-theta));
dphase = -i*(opt.gamma+opt.delta)*phase/opt.scale;

w = dfw.*phase + fw.*dphase;
w = w/sqrt(opt.scale);

%tol = 1e-6;

%phase = exp(i*(opt.gamma+opt.delta)/2*(pi-theta)).*...
%        2^((opt.gamma+opt.delta-4)/2);

%if abs(opt.delta)>tol
%  phase = phase .* abs(sin(theta/2).^(opt.delta-1));
%end
%if abs(opt.gamma)>tol
%  phase = phase .* cos(theta/2).^(opt.gamma-1);
%end
%
%w = phase .* (opt.delta*(1+exp(-i*theta)) - opt.gamma*(1+exp(i*theta)));

% Weight normalization + derivative Jacobian
% >>>>>> is this right???
%w = w/opt.scale.^2;
