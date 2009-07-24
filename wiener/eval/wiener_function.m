function[w] = wiener_function(x,k,varargin)
% [W] = WIENER_FUNCTION(X,K,{S=1, T=0, SHIFT=0, SCALE=1})
%     
%     Evaluates the Wiener function \Phi_K(x) of class (S,T), which is a direct
%     map of the generalized Fourier series (gamma=S-1,delta=T). The input X is
%     any real number, and SHIFT and SCALE dictate the affine scaling of the
%     functions. The output W has size length(X) x length(K). 

global handles;
fourier = handles.speclab.fourier;
rx = handles.speclab.wiener.maps;
opt = handles.speclab.wiener.defaults(varargin{:});

x = x(:);
theta = rx.x_to_theta(x,opt);

% The `unweighted' Wiener functions are a direct map of the generalized Fourier
% series:
fopt = struct('gamma',s-1,'delta',t);
w = fourier.eval.fseries(theta,k,fopt);
