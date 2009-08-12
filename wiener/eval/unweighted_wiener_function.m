function[w] = unweighted_wiener_function(x,k,varargin)
% [w] = unweighted_wiener_function(x,k,{s=1, t=0, shift=0, scale=1})
%     
%     Evaluates the Wiener function \Phi_k(x) of class (s,t), which is a direct
%     map of the generalized Fourier series (gamma=s-1,delta=t). The input X is
%     any real number, and shift and scale dictate the affine scaling of the
%     functions. The output w has size length(x) x length(k). 

global handles;
fourier = handles.speclab.fourier;
rx = handles.speclab.wiener.maps;
opt = handles.speclab.wiener.defaults(varargin{:});

% Force column vector
x = x(:);
theta = rx.x_to_theta(x,opt);

% The `unweighted' Wiener functions are a direct map of the generalized Fourier
% series:
fopt = struct('gamma',opt.s-1,'delta',opt.t);
w = fourier.eval.fseries(theta,k,fopt);
