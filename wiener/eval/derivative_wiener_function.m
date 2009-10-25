function[w] = derivative_wiener_function(x,k,varargin)
% [w] = derivative_wiener_function(x,k,{s=1, t=0, shift=0, scale=1})
% 
%     Evaluates the derivative of the generalized Wiener rational functions
%     \phi_k^(x) of class (s,t). These are a mapping and a weighting of the
%     canonical Fourier series. The parameters shift and scale dictate the
%     affine mapping. The output w has size length(x) x length(k).

global packages;
opt = packages.speclab.wiener.defaults(varargin{:});
wiener = packages.speclab.wiener;

x = x(:);
weight = wiener.weights.phase_shifted_sqrt_weight(x,opt);
dweight = wiener.weights.dphase_shifted_sqrt_weight(x,opt);

w = wiener.eval.unweighted_wiener_function(x,k,opt);
dw = wiener.eval.derivative_unweighted_wiener_function(x,k,opt);

N = length(x);
w = spdiags(dweight,0,N,N)*w + ...
    spdiags(weight,0,N,N)*dw;
