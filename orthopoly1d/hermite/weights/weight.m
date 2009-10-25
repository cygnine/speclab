function[w] = weight(x,varargin)
% weight -- evalutes the Hermite polynomial weight function
%
% [w] = weight(r,{mu=0, scale=1, shift=0})
%
%     Evaluates the weight function assocaited with the Hermite polynomials. The
%     weight function incorporates the affine scale and shift. (Including the
%     affine Jacobian)

global packages;
opt = packages.speclab.orthopoly1d.hermite.defaults(varargin{:});
sss = packages.speclab.common.standard_scaleshift_1d;
x = sss(x,opt);

w = (x.^2).^opt.mu.*exp(-x.^2);
w = w/opt.scale;
