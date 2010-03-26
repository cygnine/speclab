function[w] = weight(theta,varargin)
% function[w] = weight(theta,{gamma=0, delta=0, scale=1, shift=0})
%
%     Evaluates the weight function for the generalized Fourier Series. The
%     parameters gamma and delta specify which function class to consider. The
%     affine shift and scale parameters indicate the interval over which to
%     evaluate the weight function. Note that this weight function depends on
%     scale. 

persistent sss defaults weight theta_to_r
if isempty(defaults)
  from speclab.common import standard_scaleshift_1d as sss
  from speclab.fourier import defaults
  from speclab.orthopoly.jacobi.weights import weight
  from speclab.fourier.maps import theta_to_r
end

opt = defaults(varargin{:});

%theta = sss(theta,opt);
%w = ((1-cos(theta)).^opt.delta) .* ((1+cos(theta)).^opt.gamma);

r = theta_to_r(theta,opt);
w = weight(r,'alpha',opt.delta, 'beta', opt.gamma);
w = w/opt.scale;
