function[theta,w] = gauss_quadrature(N,varargin)
% [theta,wr] = gauss_quadrature(N,{gamma=0, delta=0, shift=0, scale=1})
%
%     Returns N-point the Gauss-Quadrature nodes + weights for the generalized
%     Fourier Series. The default interval is [-pi,pi].
%
%     The parameters gamma, delta determine the weight function for which this
%     quadrature rule is valid. gamma=delta=0 is the canonical Fourier set.
%     shift and scale are affine scaling parameters. 

persistent gq grq pss defaults
if isempty(pss)
  from speclab.common import physical_scaleshift_1d as pss
  from speclab.fourier import defaults
  from speclab.orthopoly.jacobi.quad import gauss_quadrature as gq
  from speclab.orthopoly.jacobi.quad import gauss_radau_quadrature as grq
end

opt = defaults(varargin{:});

if mod(N,2)==0

  [r,wr] = gq(N/2, 'alpha', opt.delta-1/2,...
                                        'beta',  opt.gamma-1/2);
  r = flipud(r); wr = flipud(wr);

  r = acos(r);
  theta = [-flipud(r); r];
  w = [flipud(wr); wr];
  theta = pss(theta,opt);

else
  
  [r,wr] = grq((N+1)/2, 'alpha', opt.delta-1/2, ...
                                                  'beta',  opt.gamma-1/2, ...
                                                  'r', 1);

  r = flipud(r); wr = flipud(wr);
  wr(1) = wr(1)*2;

  r = acos(r);
  theta = [-flipud(r); r(2:end)];
  w = [flipud(wr); wr(2:end)];
  theta = pss(theta,opt);
end
