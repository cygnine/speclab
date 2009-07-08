function[x,w] = laguerrefun_pgrq(N,alpha,xa,shift,scale)

% [x,w] = laguerrefun_pgq(N,alpha,a,shift,scale)
% Returns the N-point pi-Gauss-Radau quadrature rule for the Laguerre functions.
% The fixed node is located at x = xa;
%
% 20080623: acn

laguerre_parameters;

[x,w] = laguerrepoly_grq(N,alpha,xa,shift,scale);

%[as,bs] = laguerre_recurrence(N+1,alpha,shift,scale);
%[x,w] = opoly_gq(as,bs,N);

w = w./wlaguerre(x,alpha,shift,scale);
