function[p] = dlaguerrepolyn(x,n,alpha,shift,scale)

% [p] = dlaguerrepolyn(x,n,alpha,shift,scale)
% Evaluates the derivative of the normalized nth generalized Laguerre polynomials at the locations x. 
% The weight function is (1/scale*(t-shift))^alpha*exp(-1/scale*(t-shift))
% Is vectorized in x and n.
%
% 20080623: acn

laguerre_parameters;

N = max(n)+2;

[as,bs] = laguerre_recurrence(N+1,alpha,shift,scale);

p = eval_dopolyn(x,as,bs,n);
