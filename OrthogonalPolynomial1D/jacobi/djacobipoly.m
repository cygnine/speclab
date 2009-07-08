function[p] = djacobipoly(x,n,alpha,beta,shift,scale);

% [p] = djacobipoly(x,n,alpha,beta,shift,scale)
% Evaluates the derviatve of the monic nth generalized Jacobi polynomials at the locations x. 
% The weight function is (1-1/scale*(x-shift))^alpha*(1+1/scale*(x-shift))^beta
% Is vectorized in x and n.
%
% 20080524: acn

jacobi_parameters;

N = max(n)+1;

[as,bs] = jacobi_recurrence(N,alpha,beta,shift,scale);

p = eval_dopoly(x,as,bs,n);
