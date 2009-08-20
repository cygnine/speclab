function[g] = jacobiconst_gammatilde(ns,a,b);

% [g] = jacobiconst_gammatilde(n,a,b)
% Returns the jacobi constant \tilde{g}^(a,b)_n, which is the recurrence constant for
% (1+-r)*P^(a,b)_n, the normalized polynomials.
% Is vectorized in n.
%
% (1-r)*P_n^(a,b) = g(1)*P_(n+1)^(a-1,b) + g(2)*P_n^(a-1,b)
%
% 20080712 -- acn

ns = ns(:);
g = zeros([length(ns) 2]);

temp = 2*ns+a+b;
g(:,1) = sqrt(2*(ns+a).*(ns+a+b)./(temp.*(temp+1)));
g(:,2) = sqrt(2*(ns+1).*(ns+b+1)./((temp+1).*(temp+2)));

tol = 1e-12;
n0 = (ns==0);
if (any(n0))&&(abs(a+b)<tol);
  g(n0,1) = sqrt(2*a./(a+b+1));
end
