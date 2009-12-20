function[dw] = derivative_unweighted_wiener_function(x,k,varargin)
% [dw] = derivative_unweighted_wiener_function(x,k,{s=1, t=0, shift=0, scale=1})
%     
%     Evaluates the derivative of the Wiener function \Phi_k(x) of class (s,t),
%     which is a direct map of the generalized Fourier series
%     (gamma=s-1,delta=t). The input X is any real number, and shift and scale
%     dictate the affine scaling of the functions. The output w has size
%     length(x) x length(k). 

persistent x_to_theta dtheta_dx dfseries defaults
if isempty(defaults)
  from speclab.wiener import defaults
  from speclab.fourier.eval import dfseries
  from speclab.wiener.maps import x_to_theta
  from speclab.wiener.maps import dtheta_dx
end

opt = defaults(varargin{:});

% Force column vector
x = x(:);
N = length(x);
theta = x_to_theta(x,opt);
dtdx = dtheta_dx(x,opt);

fopt = struct('gamma',opt.s-1,'delta',opt.t);
dw = spdiags(dtdx,0,N,N)*dfseries(theta,k,fopt);
