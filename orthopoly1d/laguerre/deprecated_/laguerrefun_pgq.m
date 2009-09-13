function[x,w] = laguerrefun_pgq(N,alpha,shift,scale)

% [x,w] = laguerrefun_pgq(N,alpha,shift,scale)
% Returns the N-point pi-Gaussian quadrature rule for the Laguerre functions.
%
% 20080623: acn

laguerre_parameters;

[as,bs] = laguerre_recurrence(N+1,alpha,shift,scale);

[x,w] = opoly_gq(as,bs,N);

w = w./wlaguerre(x,alpha,shift,scale);
