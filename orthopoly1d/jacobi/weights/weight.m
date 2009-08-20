function[w] = weight(r,varargin)
% [w] = weight(r,{alpha=-1/2,beta=-1/2,shift=0,scale=1})
%
%     Evaluates the Jacobi polynomial weight function for class (alpha,beta)
%     under the affine map r -> r*scale + shift. Note that this function also
%     takes into account the Jacobian of the transformation so that if scale is
%     not 1, then the polynomials that were once 'normal', aren't anymore.

global handles;
opt = handles.speclab.orthopoly1d.jacobi.defaults(varargin{:});
sss = handles.speclab.common.standard_scaleshift_1d;
r = sss(r,opt);

w = (1-r).^opt.alpha;
w = w.*(1+r).^opt.beta;
w = w/opt.scale;
