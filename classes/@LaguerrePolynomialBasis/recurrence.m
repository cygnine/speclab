function[a,b] = recurrence(self, n)
% [a,b] = recurrence(n,d=1);
%
%     Returns the recurrence coefficients with indices n for the Laguerre
%     polynomials defined by alpha

a = zeros(size(n));
b = zeros(size(n));

a = 2*n + self.alpha + 1;

neq0 = (n==0);
b(neq0) = gamma(1 + self.alpha);
b(~neq0) = n(~neq0).*(n(~neq0) + self.alpha);
