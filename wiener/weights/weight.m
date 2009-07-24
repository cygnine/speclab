function[w] = weight(x,varargin)
% [W] = WEIGHT(X,{S=1, T=0, SCALE=1, SHIFT=0})
%
%     Evaluates the weight function for the (unweighted) Wiener rational
%     functions. The parameters S and T specify which function class to
%     consider, and the affine parameters SHIFT and SCALE determine the affine
%     map of the functions. This weight function does depend on SCALE. 

global handles;
wf = handles.speclab.fourier.weights.weight;
opt = handles.speclab.wiener.defaults(varargin{:});

theta = rx.x_to_theta(x,opt);
% Note: it just so happens that dtheta/dx = (1-cos(theta)). Therefore, we
% implicitly take default gamma=s-1 to gamma = (s-1)+1 = s.
fopt = struct('gamma',s,'delta',t);

w = wf(theta,fopt);
w = w/opt.scale;
