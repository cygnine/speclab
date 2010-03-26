function[p] = hermitefunn(x,n,mu,shift,scale)

% [p] = hermitefunn(x,n,mu,shift,scale)
% Evaluates the normalized nth generalized Hermite functions at the locations x. 
% They are orthonormal over the real line with unit weight
% Is vectorized in x and n.
%
% 20080624: acn

hermite_parameters;

N = max(n)+2;

[as,bs] = hermite_recurrence(N+1,mu,shift,scale);

p = eval_opolyn(x,as,bs,n);

X = length(x);
p = spdiags(whermite_sqrt(x,mu,shift,scale),0,X,X)*p;
