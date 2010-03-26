function[p] = jacobifun(x,n,alpha,beta,shift,scale);

% [p] = jacobifun(x,n,alpha,beta,shift,scale)
% Evaluates the monic nth generalized Jacobi functions at the locations x. 
% They are orthogonal over [-1,1] with unit weight
% Is vectorized in x and n.
%
% 20080524: acn

jacobi_parameters;

N = max(n)+1;

[as,bs] = jacobi_recurrence(N,alpha,beta,shift,scale);

p = eval_opoly(x,as,bs,n);

X = size(x,1);
w = wjacobi_sqrt(x,alpha,beta,shift,scale);

p = spdiags(w,0,X,X)*p;
