function[p] = dlaguerrefun(x,n,alpha,shift,scale)

% [p] = dlaguerrefun(x,n,alpha,shift,scale)
% Evaluates the derivative of the monic nth generalized Laguerre function at the locations x. 
% Is vectorized in x and n.
%
% 20080624: acn

laguerre_parameters;

N = max(n)+1;
x = x(:);

[as,bs] = laguerre_recurrence(N,alpha,shift,scale);

p = eval_opoly(x,as,bs,n);
dp = eval_dopoly(x,as,bs,n);
w = wlaguerre_sqrt(x,alpha,shift,scale);
dw = dwlaguerre_sqrt(x,alpha,shift,scale);
X = size(x,1);

p = spdiags(w,0,X,X)*dp + spdiags(dw,0,X,X)*p;
