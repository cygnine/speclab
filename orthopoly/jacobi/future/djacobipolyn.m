function[p] = djacobipolyn(x,n,alpha,beta,shift,scale)

% [p] = djacobipolyn(x,n,alpha,beta,shift,scale)
% Evaluates the derivative of the normalized nth generalized Jacobi polynomials at the locations x. 
% The weight function is (1-1/scale*(x-shift))^alpha*(1+1/scale*(x-shift))^beta
% Is vectorized in x and n.
%
% 20080524: acn

jacobi_parameters;

N = max(n)+2;

[as,bs] = jacobi_recurrence(N+1,alpha,beta,shift,scale);

p = eval_dopolyn(x,as,bs,n);
