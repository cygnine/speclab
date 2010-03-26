function[w] = weight(x,varargin)
% weight -- evalutes the Laguerre polynomial weight function
%
% [w] = weight(x,{alpha=0, scale=1, shift=0})
%
%     Evaluates the weight function assocaited with the Laguerre polynomials. The
%     weight function incorporates the affine scale and shift. (Including the
%     affine Jacobian)

persistent defaults sss
if isempty(defaults)
  from speclab.orthopoly.laguerre import defaults
  from speclab.common import standard_scaleshift_1d as sss
end

opt = defaults(varargin{:});
x = sss(x,opt);

w = x.^opt.alpha.*exp(-x);
w = w/opt.scale;
