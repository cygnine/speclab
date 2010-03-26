function[w] = weight(r,varargin)
% [w] = weight(r,{alpha=-1/2,beta=-1/2,shift=0,scale=1,weight_normalization=[]})
%
%     Evaluates the Jacobi polynomial weight function for class (alpha,beta)
%     under the affine map r -> r*scale + shift. When scaled, this weight
%     function absorbs the Jacobian, leaving the point evaluations of the
%     polynomials unchanged with respect to the input 'scale'.

persistent defaults sss npre recur
if isempty(defaults)
  from speclab.orthopoly.jacobi import defaults
  from speclab.orthopoly.jacobi import recurrence as recur
  from speclab.common import standard_scaleshift as sss
  from speclab.common.tensor import node_preprocessing as npre
end

opt = defaults(varargin{:});

[r,garbage] = npre(r, opt.dim);
r = sss(r,opt);

w = ones([size(r,1) 1]);
for q = 1:opt.dim 
  w = w.*(1-r(:,q)).^opt.alpha(q);
  w = w.*(1+r(:,q)).^opt.beta(q);
end

switch opt.weight_normalization
case 'probability'
  for q = 1:opt.dim
    [a,b] = recur(1,'alpha',opt.alpha(q),'beta',opt.beta(q));
    w = w/(b*opt.scale(q));
  end
otherwise
  w = w/opt.scale;
end
