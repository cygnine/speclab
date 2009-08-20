function[p] = hermitepolyn(x,n,mu,shift,scale)

% [p] = hermitepolyn(x,n,mu,shift,scale)
% Evaluates the normalized nth generalized Hermite polynomials at the locations x. 
% The weight function is \abs{(t-shift)/scale}^{2\mu} e^{-((t-shift)/scale)^2}. 
% Is vectorized in x and n.
%
% 20080524: acn

hermite_parameters;

N = max(n)+2;

[as,bs] = hermite_recurrence(N+1,mu,shift,scale);

p = eval_opolyn(x,as,bs,n);
