function[x,w] = jacobifun_pgq(N,alpha,beta,shift,scale)

% [x,w] = jacobifun_pgq(N,alpha,beta.shift,scale)
% Returns the N-point pi-Gaussian quadrature rule for the Jacobi functions with
%
% 20080623: acn

jacobi_parameters;

[as,bs] = jacobi_recurrence(N+1,alpha,beta,shift,scale);

[x,w] = opoly_gq(as,bs,N);

w = w./wjacobi(x,alpha,beta,shift,scale);
