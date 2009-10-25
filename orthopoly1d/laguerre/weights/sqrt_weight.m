function[w] = sqrt_weight(x,varargin)
% sqrt_weight -- evalutes the Laguerre polynomial weight function
%
% [w] = sqrt_weight(x,{alpha=0, scale=1, shift=0})
%
%     Evaluates the square root of the weight function assocaited with the
%     Laguerre polynomials. The weight function incorporates the affine scale
%     and shift. (Including the affine Jacobian)

global handles;
opt = handles.speclab.orthopoly1d.laguerre.defaults(varargin{:});
sss = handles.speclab.common.standard_scaleshift_1d;
x = sss(x,opt);

w = x.^(opt.alpha/2).*exp(-x/2);
w = w/sqrt(opt.scale);
