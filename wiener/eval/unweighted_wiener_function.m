function[w] = unweighted_wiener_function(x,k,varargin)
% [w] = unweighted_wiener_function(x,k,{s=1, t=0, shift=0, scale=1})
%     
%     Evaluates the Wiener function \Phi_k(x) of class (s,t), which is a direct
%     map of the generalized Fourier series (gamma=s-1,delta=t). The input X is
%     any real number, and shift and scale dictate the affine scaling of the
%     functions. The output w has size length(x) x length(k). 

persistent defaults x_to_theta fseries
if isempty(defaults)
  from speclab.wiener import defaults
  from speclab.fourier.eval import fseries
  from speclab.wiener.maps import x_to_theta
end

opt = defaults(varargin{:});

% Force column vector
x = x(:);
theta = x_to_theta(x,opt);

% The `unweighted' Wiener functions are a direct map of the generalized Fourier
% series:
fopt = struct('gamma',opt.s-1,'delta',opt.t);
w = fseries(theta,k,fopt);
