function[theta,w] = gauss_quadrature(N,varargin)
% [THETA,WR] = GAUSS_QUADRATURE(N,{GAMMA=0, DELTA=0, SHIFT=0, SCALE=1})
%
%     Returns the Gauss-Quadrature nodes + weights for the generalized Fourier
%     Series. The default interval is [0,2*pi].
%
%     The parameters GAMMA, DELTA determine the weight function for which this
%     quadrature rule is valid. GAMMA=DELTA=0 is the canonical Fourier set.
%     SHIFT and SCALE are affine scaling parameters. 

global handles;
opt = handles.speclab.fourier.defaults(varargin{:});
jac = handles.speclab.OrthogonalPolynomial1D.jacobi;
pss = handles.speclab.common.physical_scaleshift_1d;

if mod(N,2)==0

  [r,wr] = jac.quad.gauss_quadrature(N/2, 'alpha', opt.delta-1/2,...
                                        'beta',  opt.gamma-1/2);
  r = flipud(r); wr = flipud(wr);

  r = acos(r);
  theta = [-flipud(r); r];
  w = [flipud(wr); wr];
  theta = pss(theta,opt);

else
  
  [r,wr] = jac.quad.gauss_radau_quadrature((N+1)/2, 'alpha', opt.delta-1/2, ...
                                                  'beta',  opt.gamma-1/2, ...
                                                  'r', 1);

  r = flipud(r); wr = flipud(wr);
  wr(1) = wr(1)*2;

  r = acos(r);
  theta = [-flipud(r); r(2:end)];
  w = [flipud(wr); wr(2:end)];
  theta = pss(theta,opt);
end
