function[a,b] = recurrence(self, n, dim)
% [a,b] = recurrence(n,d=1);
%
%     Returns the recurrence coefficients with indices n for the Hermite
%     polynomials defined by (alpha(dim), beta(dim)). 

a = zeros(size(n));
b = zeros(size(n));

neq0 = (n==0);
nodd = mod(n,2);
b(neq0) = gamma(self.mu+1/2);
b(~neq0) = 1/2*(n(~neq0));
b(nodd) = b(nodd) + self.mu;

%a = zeros([N 1]);
%b = a;
%
%b(1) = gamma(self.mu+1/2);
%oddk = 1:2:(N-1);
%b(2:end) = 1/2*(1:(N-1));
%b(oddk+1) = b(oddk+1)+self.mu;
