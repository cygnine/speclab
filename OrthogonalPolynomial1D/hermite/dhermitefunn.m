function[p] = dhermitefunn(x,n,mu,shift,scale);

% [p] = dhermitefunn(x,n,mu,shift,scale)
% Evaluates the derviatve of the normalized nth generalized Hermite function at the locations x. 
% Is vectorized in x and n.
%
% 20080624: acn

hermite_parameters;

N = max(n)+2;

[as,bs] = hermite_recurrence(N,mu,shift,scale);

p = eval_opolyn(x,as,bs,n);
dp = eval_dopolyn(x,as,bs,n);

X = length(x);
w = whermite_sqrt(x,mu,shift,scale);
dw = dwhermite_sqrt(x,mu,shift,scale);

p = spdiags(w,0,X,X)*dp + spdiags(dw,0,X,X)*p;
