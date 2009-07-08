function[p] = laguerrefunn(x,n,alpha,shift,scale)

% [p] = laguerrefunn(x,n,alpha,shift,scale)
% Evaluates the normalized nth generalized Laguerre functions at the locations x. 
% The functions are orthogonal under unit weight over
% [shift,\infty)
% Is vectorized in x and n.
%
% 20080623: acn

laguerre_parameters;

x = x(:);
N = max(n)+2;

[as,bs] = laguerre_recurrence(N+1,alpha,shift,scale);

p = eval_opolyn(x,as,bs,n);
w = wlaguerre_sqrt(x,alpha,shift,scale);
X = size(x,1);

p = spdiags(w,0,X,X)*p;
