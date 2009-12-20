function[w] = weight(x,varargin)
% [w] = weight(x,{s=1, t=0, scale=1, shift=0})
%
%     Evaluates the weight function for the (unweighted) Wiener rational
%     functions. The parameters S and T specify which function class to
%     consider, and the affine parameters shift and scale determine the affine
%     map of the functions. This weight function does depend on scale. 

persistent defaults x_to_theta wf
if isempty(defaults)
  from speclab.wiener import defaults
  from speclab.wiener.maps import x_to_theta
  from speclab.fourier.weights import weight as wf
end

opt = defaults(varargin{:});

theta = x_to_theta(x,opt);
% Note: it just so happens that dtheta/dx = (1-cos(theta)). Therefore, we
% implicitly take default gamma=s-1 to gamma = (s-1)+1 = s.
fopt = struct('gamma',opt.s,'delta',opt.t);

w = wf(theta,fopt);
w = w/opt.scale;
