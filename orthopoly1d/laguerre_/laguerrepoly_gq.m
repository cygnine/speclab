function[x,w] = laguerrepoly_gq(N,alpha,shift,scale)

% [x,w] = laguerrepoly_gq(N,alpha,shift,scale)
% Returns the N-point Gaussian quadrature rule for the Laguerre polynomials with
% weight function (1/scale*(t-shift))^alpha*exp(-1/scale*(t-shift))
%
% 20080623: acn

laguerre_parameters;

[as,bs] = laguerre_recurrence(N+1,alpha,shift,scale);

[x,w] = opoly_gq(as,bs,N);
