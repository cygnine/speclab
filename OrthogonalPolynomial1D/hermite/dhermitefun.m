function[p] = dhermitefun(x,n,mu,shift,scale);

% [p] = dhermitefun(x,n,mu,shift,scale)
% Evaluates the derviatve of the monic nth generalized Hermite function at the locations x. 
% Is vectorized in x and n.
%
% 20080624: acn

hermite_parameters;

N = max(n)+1;

[as,bs] = hermite_recurrence(N,mu,shift,scale);

p = eval_opoly(x,as,bs,n);
dp = eval_dopoly(x,as,bs,n);

X = length(x);
w = whermite_sqrt(x,mu,shift,scale);
dw = dwhermite_sqrt(x,mu,shift,scale);

p = spdiags(w,0,X,X)*dp + spdiags(dw,0,X,X)*p;
