function[h] = jacobiconst_etatilde(ns,a,b);

% [h] = jacobiconst_etatilde(n,a,b)
% Returns the jacobi constant eta^(a,b)_n, which is the recurrence constant for
% P^(a,b)_n, the normalized polynomials.
% Is vectorized in n.
%
% P^(a,b)_n = h(1)*P^(a+1,b+1)_n + h(2)*P^(a+1,b+1)_(n-1) +
% h(3)*P^(a+1,b+1)_(n-2).
%
% 20080714 -- acn

ns = ns(:);
h = zeros([length(ns) 3]);

denom = (2*ns+a+b);
h(:,1) = sqrt(4*(ns+a+1).*(ns+b+1).*(ns+a+b+1).*(ns+a+b+2)./...
((denom+1).*(denom+2).^2.*(denom+3)));
h(:,2) = (2*(a-b)*sqrt(ns.*(ns+a+b+1)))./(denom.*(denom+2));
h(:,3) = -sqrt(4*ns.*(ns-1).*(ns+a).*(ns+b)./((denom-1).*denom.^2.*(denom+1)));
