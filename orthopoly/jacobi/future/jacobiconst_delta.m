function[d] = jacobiconst_delta(ns,a,b);

% [d] = jacobiconst_delta(n,a,b)
% Returns the jacobi constant d^(a,b)_n, which is the recurrence constant for
% P^(a,b)_n, the monic polynomials.
% Is vectorized in n.
%
% P^(a,b)_n = d(n,a,b)(1)*P_n^(a+1,b) - d(n,a,b)(2)*P^(a+1,b)_(n-1)
% P^(a,b)_n = d(n,b,a)(1)*P_n^(a,b+1) + d(n,b,a)(2)*P^(a,b+1)_(n-1)
%
% 20080712 -- acn

ns = ns(:);
d = zeros([length(ns) 2]);

d(:,1) = 1;
d(:,2) = 2*ns.*(ns+b)./((2*ns+a+b).*(2*ns+a+b+1));
