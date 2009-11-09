function[w] = weight(r,varargin)
% [w] = weight(r,{alpha=-1/2,beta=-1/2,shift=0,scale=1,weight_normalization=[]})
%
%     Evaluates the Jacobi polynomial weight function for class (alpha,beta)
%     under the affine map r -> r*scale + shift. Note that this function also
%     takes into account the Jacobian of the transformation so that if scale is
%     not 1, then the polynomials that were once 'normal', aren't anymore.

from speclab.orthopoly1d.jacobi defaults
sss = from_as('speclab.common', 'standard_scaleshift_1d');
opt = defaults(varargin{:});
r = sss(r,opt);

w = (1-r).^opt.alpha;
w = w.*(1+r).^opt.beta;

switch opt.weight_normalization
case 'probability'
  recur = from_as('speclab.orthopoly1d.jacobi.coefficients', 'recurrence');
  [a,b] = recur(1,'alpha',opt.alpha,'beta',opt.beta);
  w = w/(b*opt.scale);
otherwise
  w = w/opt.scale;
end
