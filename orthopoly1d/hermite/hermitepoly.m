function[p] = hermitepoly(x,n,mu,shift,scale)

% [p] = hermitepoly(x,n,mu,shift,scale)
% Evaluates the monic nth generalized Hermite polynomials at the locations x. 
% The weight function is \abs{(t-shift)/scale}^{2\mu} e^{-((t-shift)/scale)^2}. 
% Is vectorized in x and n.
%
% 20080524: acn

hermite_parameters;

N = max(n)+1;

[as,bs] = hermite_recurrence(N,mu,shift,scale);

p = eval_opoly(x,as,bs,n);
