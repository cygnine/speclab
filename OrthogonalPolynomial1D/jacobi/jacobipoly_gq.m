function[x,w] = jacobipoly_gq(N,alpha,beta,shift,scale)

% [x,w] = jacobipoly_gq(N,alpha,beta.shift,scale)
% Returns the N-point Gaussian quadrature rule for the Jacobi polynomials with
% The weight function is (1-1/scale*(x-shift))^alpha*(1+1/scale*(x-shift))^beta
%
% 20080623: acn

jacobi_parameters;

tol = 1e-8;
if (abs(alpha+1/2)<tol)&&(abs(beta+1/2)<tol);
  w = pi/N*ones([N 1]);
  temp = linspace(pi,0,N+1)';
  x = cos(temp(1:N)-pi/(2*N));
else
  [as,bs] = jacobi_recurrence(N+1,alpha,beta,shift,scale);

  [x,w] = opoly_gq(as,bs,N);
end
