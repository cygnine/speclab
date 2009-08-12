function[w] = dweight(r,varargin)
% [w] = dweight(r,{alpha=-1/2,beta=-1/2,shift=0,scale=1})
%
%     Evaluates the derivative of the Jacobi polynomial weight function for
%     class (alpha,beta) under the affine map r -> r*scale + shift. Note that
%     this function also takes into account the Jacobian of the transformation
%     so that if scale is not 1, then the polynomials that were once 'normal',
%     aren't anymore.

global handles;
opt = handles.speclab.OrthogonalPolynomial1D.jacobi.defaults(varargin{:});
sss = handles.speclab.common.standard_scaleshift_1d;
r = sss(r,opt);

wa = (1-r).^opt.alpha;
if abs(opt.alpha)>1e-8
  dwa = opt.alpha*(1-r).^(opt.alpha-1);
else
  dwa = 0;
end

wb = (1+r).^opt.beta;
if abs(opt.beta)>1e-8
  dwb = opt.beta*(1+r).^(opt.beta-1);
else
  dwb = 0;
end

w = dwa.*wb + wa.*dwb;
w = w/opt.scale^2;
