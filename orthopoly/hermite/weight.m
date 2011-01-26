function[w] = weight(x,varargin)
% weight -- evalutes the Hermite polynomial weight function
%
% [w] = weight(r,{mu=0, scale=1, shift=0})
%
%     Evaluates the weight function assocaited with the Hermite polynomials. The
%     weight function incorporates the affine scale and shift. (Including the
%     affine Jacobian)

persistent defaults sss recurrence npre
if isempty(defaults)
  from speclab.orthopoly.hermite import defaults recurrence
  from speclab.common import standard_scaleshift as sss
  from speclab.common.tensor import node_preprocessing as npre
end

opt = defaults(varargin{:});
[x, garbage] = npre(x, opt.dim);
x = sss(x,opt);

w = ones([size(x, 1) 1]);
for q = 1:opt.dim
  w = w.*(x(:,q).^2).^opt.mu(q).*exp(-x(:,q).^2);
end

switch opt.weight_normalization
case 'probability'
  for q = 1:opt.dim
    [a,b] = recurrence(1,'mu',opt.mu(q));
    w = w/(b*opt.scale(q));
  end
otherwise
  w = w/prod(opt.scale);
end
