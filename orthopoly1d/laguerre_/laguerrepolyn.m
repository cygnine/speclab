function[p] = laguerrepolyn(x,n,alpha,shift,scale)

% [p] = laguerrepolyn(x,n,alpha,shift,scale)
% Evaluates the normalized nth generalized Laguerre polynomials at the locations x. 
% The weight function is (1/scale*(t-shift))^alpha*exp(-1/scale*(t-shift))
% Is vectorized in x and n.
%
% 20080623: acn

laguerre_parameters;

N = max(n)+2;

[as,bs] = laguerre_recurrence(N+1,alpha,shift,scale);

p = eval_opolyn(x,as,bs,n);
