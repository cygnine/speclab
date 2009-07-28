function[w] = weight(r,varargin)
% [W] = WEIGHT(R,{ALPHA=-1/2,BETA=-1/2,SHIFT=0,SCALE=1})
%
%     Evaluates the Jacobi polynomial weight function for class (ALPHA,BETA)
%     under the affine map r -> r*SCALE + SHIFT. Note that this function also
%     takes into account the Jacobian of the transformation so that if SCALE is
%     not 1, then the polynomials that were once 'normal', aren't anymore.

global handles;
opt = handles.speclab.OrthogonalPolynomial1D.jacobi.defaults(varargin{:});

w = (1-r).^opt.alpha;
w = w.*(1+r).^opt.beta;
w = w/opt.scale;
