function[w] = weight(theta,varargin)
% function[w] = weight(theta,{gamma=0, delta=0, scale=1, shift=0})
%
%     Evaluates the weight function for the generalized Fourier Series. The
%     parameters gamma and delta specify which function class to consider. The
%     affine shift and scale parameters indicate the interval over which to
%     evaluate the weight function. Note that this weight function depends on
%     scale. 

global handles;
sss = handles.speclab.common.standard_scaleshift_1d;
opt = handles.speclab.fourier.defaults(varargin{:});
jac = handles.speclab.OrthogonalPolynomial1D.jacobi;
xr = handles.speclab.fourier.maps;

%theta = sss(theta,opt);
%w = ((1-cos(theta)).^opt.delta) .* ((1+cos(theta)).^opt.gamma);

r = xr.theta_to_r(theta,opt);
w = jac.weights.weight(r,'alpha',opt.delta, 'beta', opt.gamma);
w = w/opt.scale;
