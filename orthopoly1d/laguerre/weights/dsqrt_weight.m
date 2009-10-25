function[dw] = dsqrt_weight(x,varargin)
% dsqrt_weight -- evalutes the Laguerre polynomial weight function
%
% [dw] = dsqrt_weight(x,{alpha=0, scale=1, shift=0})
%
%     Evaluates the derivative of the square root of the weight function
%     assocaited with the Laguerre polynomials. The weight function incorporates
%     the affine scale and shift. (Including the affine Jacobian)

global handles;
opt = handles.speclab.orthopoly1d.laguerre.defaults(varargin{:});
sss = handles.speclab.common.standard_scaleshift1d;
x = sss(x,opt);

if opt.alpha==0
  dw = -1/2*exp(-x/2);
else
  dw = (opt.alpha/2.*x.^(opt.alpha/2-1) - 1/2*x.^(opt.alpha/2)).*exp(-x/2);
end

dw = dw/(opt.scale^(3/2));
