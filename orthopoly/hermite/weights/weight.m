function[w] = weight(x,varargin)
% weight -- evalutes the Hermite polynomial weight function
%
% [w] = weight(r,{mu=0, scale=1, shift=0})
%
%     Evaluates the weight function assocaited with the Hermite polynomials. The
%     weight function incorporates the affine scale and shift. (Including the
%     affine Jacobian)

defaults = from_as('speclab.orthopoly.hermite', 'defaults');
sss = from_as('speclab.common', 'standard_scaleshift_1d');
opt = defaults(varargin{:});
x = sss(x,opt);

w = (x.^2).^opt.mu.*exp(-x.^2);

switch opt.weight_normalization
case 'probability'
  recur = from_as('speclab.orthopoly.hermite.coefficients', 'recurrence');
  [a,b] = recur(1,'mu',opt.mu);
  w = w/(b*opt.scale);
otherwise
  w = w/opt.scale;
end
