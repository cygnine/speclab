function[p] = hermitefun(x,n,mu,shift,scale)

% [p] = hermitefun(x,n,mu,shift,scale)
% Evaluates the monic nth generalized Hermite functions at the locations x. 
% They are orthogonal over the real line with unit weight.
% Is vectorized in x and n.
%
% 20080624: acn

hermite_parameters;

N = max(n)+1;

[as,bs] = hermite_recurrence(N,mu,shift,scale);

p = eval_opoly(x,as,bs,n);

X = length(x);
p = spdiags(whermite_sqrt(x,mu,shift,scale),0,X,X)*p;
