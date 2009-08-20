function[h] = jacobiconst_h(ns,a,b);

% [h] = jacobiconst_h(n,a,b)
% Returns the jacobi constant h^(a,b)_n, which is the square of the 2-norm of
% the monic polynomial P^(a,b)_n. 
% Is vectorized in n.
%
% 20080712 -- acn

%h = 2^(a+b+1)./(2*n+a+b+1)*gamma(n+a+1)*gamma(n+b+1)./...
%  (gamma(n+1)*gamma(n+a+b+1));

%[as,bs] = jacobi_recurrence(a,b,max(ns(:))+1);

[as,bs] = jacobi_recurrence(max(ns(:))+1,a,b);

hs = cumprod(bs);

h = hs(ns+1);

% Actual formula (numerically accurate, but overflows quickly->NaN):
%h(:,2) = 2.^(2*ns+a+b+1).*factorial(ns).*gamma(ns+a+b+1).*gamma(ns+a+1).*...
%gamma(ns+b+1)./((2*ns+a+b+1).*gamma(2*ns+a+b+1).^2);
