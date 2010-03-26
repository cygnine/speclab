function[w] = dweight(x,varargin)
% weight -- evalutes derivative of the Hermite polynomial weight function
%
% [w] = dweight(r,{mu=0, scale=1, shift=0})
%
%     Evaluates the derivative of the weight function assocaited with the
%     Hermite polynomials. The weight function incorporates the affine scale and
%     shift. (Including the affine Jacobian)

persistent defaults sss
if isempty(defaults)
  from speclab.orthopoly.hermite import defaults
  from speclab.common import standard_scaleshift_1d as sss
end

opt = defaults(varargin{:});
x = sss(x,opt);

if mu==0
  w = -2*x.*exp(-x.^2);
else
  w = 2*abs(x).^(2*mu-1).*exp(-x.^2).*(sign(x)*mu - x.*abs(x));
end
w = w/opt.scale.^2;
