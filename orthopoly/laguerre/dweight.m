function[dw] = dweight(x,varargin)
% dweight -- evalutes the Laguerre polynomial weight function
%
% [dw] = dweight(x,{alpha=0, scale=1, shift=0})
%
%     Evaluates the derivative of the weight function assocaited with the
%     Laguerre polynomials. The weight function incorporates the affine scale
%     and shift. (Including the affine Jacobian)

persistent defaults sss
if isempty(defaults)
  from speclab.orthopoly.laguerre import defaults
  from speclab.common import standard_scaleshift_1d as sss
end

opt = defaults(varargin{:});
x = sss(x,opt);

if opt.alpha==0
  dw = -exp(-x);
else
  dw = (opt.alpha.*x.^(opt.alpha-1) - x.^opt.alpha).*exp(-x);
end

dw = dw/opt.scale.^2;
