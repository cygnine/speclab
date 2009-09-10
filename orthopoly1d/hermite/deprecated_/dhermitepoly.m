function[p] = dhermitepoly(x,n,mu,shift,scale);

% [p] = dhermitepoly(x,n,mu,shift,scale)
% Evaluates the derviatve of the monic nth generalized Hermite polynomials at the locations x. 
% The weight function is \abs{(t-shift)/scale}^{2\mu} e^{-((t-shift)/scale)^2}. 
% Is vectorized in x and n.
%
% 20080624: acn

hermite_parameters;

N = max(n)+1;

[as,bs] = hermite_recurrence(N,mu,shift,scale);

p = eval_dopoly(x,as,bs,n);
