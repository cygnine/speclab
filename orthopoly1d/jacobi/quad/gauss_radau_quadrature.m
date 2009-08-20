function[x,w] = gauss_radau_quadrature(N,varargin)
% [X,W] = GAUSS_RADAU_QUADRATURE(N,{ALPHA=-1/2,BETA=-1/2,SHIFT=0,SCALE=1,R=1})
%
%     Returns the N-point Gauss-Radau quadrature rule for the Jacobi polynomials.
%     The weight function is
%     (1-1/SCALE*(x-SHIFT))^ALPHA*(1+1/SCALE*(x-SHIFT))^BETA. The Radau point is
%     located at x=R.

global handles;
opoly = handles.speclab.orthopoly1d;
jac = opoly.jacobi;
pss = handles.speclab.common.physical_scaleshift_1d;
sss = handles.speclab.common.standard_scaleshift_1d;

opt = jac.defaults(varargin{:});
[alpha,beta,scale,shift,r] = ...
  deal(opt.alpha,opt.beta,opt.scale,opt.shift,opt.r);
r = sss(r,opt);

tol = 1e-8;
if (abs(alpha+1/2)<tol)&&(abs(beta+1/2)<tol)&&(abs(abs(r)-1)<tol);
  w = pi/N*ones([N 1]); % GQ
  w = 2*pi/(2*N-1)*ones([N 1]);
  w(1) = pi/(2*N-1);
  temp = (0:(N-1)).'; linspace(pi,0,N+1)';
  x = -cos(temp*2*pi/(2*N-1));
  if r>0;
    x = flipud(-x);
    w = flipud(w);
  end
else
  [a,b] = jac.coefficients.recurrence(N,opt);

  [x,w] = opoly.gauss_radau_quadrature(a,b,N,r);
end

x = pss(x,opt);
w = w*scale;
