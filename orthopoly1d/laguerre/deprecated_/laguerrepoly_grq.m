function[x,w] = laguerrepoly_grq(N,alpha,xa,shift,scale)

% [x,w] = laguerrepoly_grq(N,alpha,shift,scale)
% Returns the N-point Gauss-Radau quadrature rule for the Laguerre polynomials with
% the weight function (1-1/scale*(x-shift))^alpha*(1+1/scale*(x-shift))^beta
% The Radau point is at x = a
%
% 20080623: acn

laguerre_parameters;

[as,bs] = laguerre_recurrence(N,alpha,shift,scale);

[x,w] = opoly_grq(as,bs,N,xa);
