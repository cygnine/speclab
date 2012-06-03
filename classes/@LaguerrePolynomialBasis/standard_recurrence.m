function[a,b] = standard_recurrence(self, n, dim)
% [a,b] = standard_recurrence(n)
%
%     Returns the recurrence coefficients with indices n for the Hermite
%     polynomials defined by self.mu.

a = zeros(size(n));
b = zeros(size(n));

a = 2*n + self.alpha + 1;

neq0 = (n==0);
b(neq0) = gamma(1 + self.alpha);
b(~neq0) = n(~neq0).*(n(~neq0) + self.alpha);
