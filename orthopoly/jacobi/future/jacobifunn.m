function[p] = jacobifunn(x,n,alpha,beta,shift,scale)

% [p] = jacobifunn(x,n,alpha,beta,shift,scale)
% Evaluates the normalized nth generalized Jacobi functions at the locations x. 
% They are orthogonal over [-1,1] with unit weight
% Is vectorized in x and n.
%
% 20080624: acn

jacobi_parameters;

N = max(n)+2;

[as,bs] = jacobi_recurrence(N+1,alpha,beta,shift,scale);

p = eval_opolyn(x,as,bs,n);

X = size(x,1);
w = wjacobi_sqrt(x,alpha,beta,shift,scale);

p = spdiags(w,0,X,X)*p;
