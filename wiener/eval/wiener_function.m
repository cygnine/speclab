function[w] = wiener_function(x,k,varargin)
% [w] = wiener_function(x,k,{s=1, t=0, shift=0, scale=1})
% 
%     Evaluates the generalized Wiener rational functions \phi_k^(x) of class
%     (s,t). These are a mapping and a weighting of the canonical Fourier
%     series. The parameters shift and scale dictate the affine mapping. The
%     output w has size length(x) x length(k).

global handles;
opt = handles.speclab.wiener.defaults(varargin{:});

wiener = handles.speclab.wiener;
weight = wiener.weights.phase_shifted_sqrt_weight;

w = wiener.eval.unweighted_wiener_function(x,k,opt);
x = x(:);

N = length(x);
w = spdiags(weight(x,opt),0,N,N)*w;
