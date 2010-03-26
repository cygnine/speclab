function[h] = jacobiconst_eta(ns,a,b);

% [h] = jacobiconst_eta(n,a,b)
% Returns the jacobi constant eta^(a,b)_n, which is the recurrence constant for
% P^(a,b)_n, the monic polynomials.
% Is vectorized in n.
%
% P^(a,b)_n = h(1)*P^(a+1,b+1)_n + h(2)*P^(a+1,b+1)_(n-1) +
% h(3)*P^(a+1,b+1)_(n-2).
%
% 20080712 -- acn

ns = ns(:);
h = zeros([length(ns) 3]);

denom = (2*ns+a+b);
h(:,1) = 1;
h(:,2) = (2*ns*(a-b))./(denom.*(denom+2));
h(:,3) = -4*ns.*(ns-1).*(ns+a).*(ns+b)./((denom-1).*denom.^2.*(denom+1));
