function[w] = weighted_wiener_function(x,k,varargin)
% [W] = WEIGHTED_WIENER_FUNCTION(X,K,{S=1, T=0, SHIFT=0, SCALE=1})
% 
%     Evaluates the generalized Wiener rational functions \phi_k^(x) of class
%     (S,T). These are a mapping and a weighting of the canonical Fourier
%     Series. The parameters SHIFT and SCALE dictate the affine mapping. The
%     output W has size length(X) x length(K).

global handles;
opt = handles.speclab.wiener.defaults(varargin{:});

wiener = handles.speclab.wiener;
weight = wiener.weights.phase_shifted_sqrt_weight;

w = wiener.eval.wiener_function(x,k,opt);

N = length(x);
w = spdiags(weight(x,opt),0,N,N)*w;
