function[w] = wiener_function(x,k,varargin)
% [w] = wiener_function(x,k,{s=1, t=0, shift=0, scale=1})
% 
%     Evaluates the generalized Wiener rational functions \phi_k^(x) of class
%     (s,t). These are a mapping and a weighting of the canonical Fourier
%     series. The parameters shift and scale dictate the affine mapping. The
%     output w has size length(x) x length(k).

persistent defaults weight unweighted_fcn
if isempty(defaults)
  from speclab.wiener import defaults
  from speclab.wiener.weights import phase_shifted_sqrt_weight as weight
  from speclab.wiener.eval import unweighted_wiener_function as unweighted_fcn
end

opt = defaults(varargin{:});

w = unweighted_fcn(x,k,opt);
x = x(:);

N = length(x);
w = spdiags(weight(x,opt),0,N,N)*w;
