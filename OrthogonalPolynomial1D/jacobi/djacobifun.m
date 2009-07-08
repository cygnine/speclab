function[p] = djacobifun(x,n,alpha,beta,shift,scale);

% [p] = djacobifun(x,n,alpha,beta,shift,scale)
% Evaluates the derviatve of the monic nth generalized Jacobi functions at the locations x. 
% Is vectorized in x and n.
%
% 20080524: acn

jacobi_parameters;

N = max(n)+1;

[as,bs] = jacobi_recurrence(N,alpha,beta,shift,scale);

p = eval_opoly(x,as,bs,n);
dp = eval_dopoly(x,as,bs,n);
w = wjacobi_sqrt(x,alpha,beta,shift,scale);
dw = dwjacobi_sqrt(x,alpha,beta,shift,scale);
X = size(x,1);

p = spdiags(w,0,X,X)*dp + spdiags(dw,0,X,X)*p;
