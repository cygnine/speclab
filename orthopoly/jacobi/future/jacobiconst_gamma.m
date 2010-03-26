function[g] = jacobiconst_gamma(ns,a,b);

% [g] = jacobiconst_gamma(n,a,b)
% Returns the jacobi constant g^(a,b)_n, which is the recurrence constant for
% (1+-r)*P^(a,b)_n, the monic polynomials.
% Is vectorized in n.
%
% 20080712 -- acn

g = 2*(ns+a).*(ns+a+b)./((2*ns+a+b).*(2*ns+a+b+1));
