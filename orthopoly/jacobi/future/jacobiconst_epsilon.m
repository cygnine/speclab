function[e] = jacobiconst_epsilon(ns,a,b);

% [g] = jacobiconst_epsilon(n,a,b)
% Returns the jacobi constant e^(a,b)_n, which is the recurrence constant for
% (1-r^2)*P^(a,b)_n, the monic polynomials.
% Is vectorized in n.
%
% (1-r^2)*P^(a,b)_n = e(1)*P^(a-1,b-1)_n + e(2)*P^(a-1,b-1)_(n+1) +
% e(3)*P^(a-1,b-1)_(n+2).
%
% 20080712 -- acn

ns = ns(:);
e = zeros([length(ns) 3]);

denom = (2*ns+a+b);
e(:,1) = 4*(ns+a).*(ns+b).*(ns+a+b-1).*(ns+a+b)./((denom-1).*denom.^2.*(denom+1));
e(:,2) = 2*(a-b).*(ns+a+b)./(denom.*(denom+2));
e(:,3) = -1;
