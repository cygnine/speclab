function[p] = laguerrepoly(x,n,alpha,shift,scale)

% [p] = laguerrepoly(x,n,alpha,shift,scale)
% Evaluates the monic nth generalized Laguerre polynomials at the locations x. 
% The weight function is (1/scale*(t-shift))^alpha*exp(-1/scale*(t-shift)), the interval is
% [shift,\infty)
% Is vectorized in x and n.
%
% 20080623: acn

laguerre_parameters;

N = max(n)+1;

[as,bs] = laguerre_recurrence(N,alpha,shift,scale);

p = eval_opoly(x,as,bs,n);
