function[p] = djacobifunn(x,n,alpha,beta,shift,scale)

% [p] = djacobifunn(x,n,alpha,beta,shift,scale)
% Evaluates the derivative of the normalized nth generalized Jacobi functions at the locations x. 
% Is vectorized in x and n.
%
% 20080524: acn

jacobi_parameters;

N = max(n)+2;

[as,bs] = jacobi_recurrence(N+1,alpha,beta,shift,scale);

p = eval_opolyn(x,as,bs,n);
dp = eval_dopolyn(x,as,bs,n);
w = wjacobi_sqrt(x,alpha,beta,shift,scale);
dw = dwjacobi_sqrt(x,alpha,beta,shift,scale);
X = size(x,1);

p = spdiags(w,0,X,X)*dp + spdiags(dw,0,X,X)*p;
