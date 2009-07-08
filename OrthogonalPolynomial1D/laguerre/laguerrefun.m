function[p] = laguerrefun(x,n,alpha,shift,scale)

% [p] = laguerrefun(x,n,alpha,shift,scale)
% Evaluates the monic nth generalized Laguerre function at the locations x. 
% The functions are orthogonal under unit weight over
% [shift,\infty)
% Is vectorized in x and n.
%
% 20080624: acn

laguerre_parameters;

x = x(:);
N = max(n)+1;

[as,bs] = laguerre_recurrence(N,alpha,shift,scale);

p = eval_opoly(x,as,bs,n);
w = wlaguerre_sqrt(x,alpha,shift,scale);
X = size(x,1);

p = spdiags(w,0,X,X)*p;
