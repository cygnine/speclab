function[dw] = derivative_unweighted_wiener_function(x,k,varargin)
% [dw] = derivative_unweighted_wiener_function(x,k,{s=1, t=0, shift=0, scale=1})
%     
%     Evaluates the derivative of the Wiener function \Phi_k(x) of class (s,t),
%     which is a direct map of the generalized Fourier series
%     (gamma=s-1,delta=t). The input X is any real number, and shift and scale
%     dictate the affine scaling of the functions. The output w has size
%     length(x) x length(k). 

global handles;
fourier = handles.speclab.fourier;
rx = handles.speclab.wiener.maps;
opt = handles.speclab.wiener.defaults(varargin{:});

% Force column vector
x = x(:);
N = length(x);
theta = rx.x_to_theta(x,opt);
dtdx = rx.dtheta_dx(x,opt);

fopt = struct('gamma',opt.s-1,'delta',opt.t);
dw = spdiags(dtdx,0,N,N)*fourier.eval.dfseries(theta,k,fopt);
