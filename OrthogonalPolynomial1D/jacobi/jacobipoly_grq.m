function[x,w] = jacobipoly_grq(N,alpha,beta,a,shift,scale)

% [x,w] = jacobipoly_grq(N,alpha,beta,a,shift,scale)
% Returns the N-point Gauss-Radau quadrature rule for the Jacobi polynomials with
% the weight function (1-1/scale*(x-shift))^alpha*(1+1/scale*(x-shift))^beta
% The Radau point is at x = a
%
% 20080623: acn

jacobi_parameters;

tol = 1e-8;

if (abs(alpha+1/2)<tol)&&(abs(beta+1/2)<tol)&&(abs(abs(a)-1)<tol);
  w = pi/N*ones([N 1]); % GQ
  w = 2*pi/(2*N-1)*ones([N 1]);
  w(1) = pi/(2*N-1);
  temp = (0:(N-1)).'; linspace(pi,0,N+1)';
  x = -cos(temp*2*pi/(2*N-1));
  if a>0;
    x = flipud(-x);
    w = flipud(w);
  end
else
  [as,bs] = jacobi_recurrence(N,alpha,beta,shift,scale);

  [x,w] = opoly_grq(as,bs,N,a);
end
