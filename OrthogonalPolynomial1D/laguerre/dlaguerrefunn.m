function[p] = dlaguerrefunn(x,n,alpha,shift,scale)

% [p] = dlaguerrefunn(x,n,alpha,shift,scale)
% Evaluates the derivative of the normalized nth generalized Laguerre functions at the locations x. 
% Is vectorized in x and n.
%
% 20080624: acn

laguerre_parameters;

N = max(n)+2;

[as,bs] = laguerre_recurrence(N+1,alpha,shift,scale);

p = eval_opolyn(x,as,bs,n);
dp = eval_dopolyn(x,as,bs,n);
w = wlaguerre_sqrt(x,alpha,shift,scale);
dw = dwlaguerre_sqrt(x,alpha,shift,scale);
X = size(x,1);

p = spdiags(w,0,X,X)*dp + spdiags(dw,0,X,X)*p;
