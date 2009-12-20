function[w] = derivative_wiener_function(x,k,varargin)
% [w] = derivative_wiener_function(x,k,{s=1, t=0, shift=0, scale=1})
% 
%     Evaluates the derivative of the generalized Wiener rational functions
%     \phi_k^(x) of class (s,t). These are a mapping and a weighting of the
%     canonical Fourier series. The parameters shift and scale dictate the
%     affine mapping. The output w has size length(x) x length(k).

persistent defaults sqrt_weight dsqrt_weight unweighted_fun dunweighted_fun
if isempty(defaults)
  from speclab.wiener import defaults
  from speclab.wiener.weights import phase_shifted_sqrt_weight as sqrt_weight
  from speclab.wiener.weights import dphase_shifted_sqrt_weight as dsqrt_weight
  from speclab.wiener.eval import unweighted_wiener_function as unweighted_fun
  from speclab.wiener.eval import derivative_unweighted_wiener_function as dunweighted_fun
end

opt = defaults(varargin{:});

x = x(:);
weight = sqrt_weight(x,opt);
dweight = dsqrt_weight(x,opt);

w = unweighted_fun(x,k,opt);
dw = dunweighted_fun(x,k,opt);

N = length(x);
w = spdiags(dweight,0,N,N)*w + ...
    spdiags(weight,0,N,N)*dw;
