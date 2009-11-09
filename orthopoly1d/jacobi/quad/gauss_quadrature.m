function[x,w] = gauss_quadrature(N,varargin)
% [x,w] = gauss_quadrature(N,{alpha=-1/2,beta=-1/2,shift=0,scale=1})
%
%     Returns the N-point Gaussian quadrature rule for the Jacobi polynomials. 
%     The weight function is given by speclab.orthopoly1d.jacobi.weights.weight.

jac = from_as('speclab.orthopoly1d', 'jacobi');
pss = from_as('speclab.common', 'physical_scaleshift_1d');
gq = from_as('speclab.orthopoly1d', 'gauss_quadrature');

opt = jac.defaults(varargin{:});
[alpha,beta,scale,shift] = deal(opt.alpha,opt.beta,opt.scale,opt.shift);

tol = 1e-8;
if (abs(alpha+1/2)<tol)&&(abs(beta+1/2)<tol);
  w = pi/N*ones([N 1]);
  temp = linspace(pi,0,N+1)';
  x = cos(temp(1:N)-pi/(2*N));
else
  [a,b] = jac.coefficients.recurrence(N+1,opt);

  [x,w] = gq(a,b,N);
end

switch opt.weight_normalization
case 'probability'
  w = w/b(1);
end

x = pss(x,opt);
