function[x,w] = jacobipoly_glq(N,alpha,beta,shift,scale)

% [x,w] = jacobipoly_grq(N,alpha,beta,shift,scale)
% Returns the N-point Gauss-Lobatto quadrature rule for the Jacobi polynomials with
% the weight function (1-1/scale*(x-shift))^alpha*(1+1/scale*(x-shift))^beta
% The Lobatto points are at x = shift-scale and x = shift+scale
%
% 20080623: acn

jacobi_parameters;

[as,bs] = jacobi_recurrence(N,alpha,beta,shift,scale);

[x,w] = opoly_glq(as,bs,N,shift-scale,shift+scale);
