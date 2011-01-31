function[a,b] = recurrence(n)
% [a,b] = recurrence(n)
%
%     Returns the recurrence coefficients with indices n for the
%     Meixner-Pollaczek polynomials.

a = zeros(size(n));
b = zeros(size(n));

a = -(n+self.lambda)/tan(self.phi);
b = n.*(n + 2*self.lambda - 1)/(4*sin(self.phi)^2);
b(n==0) = gamma(2*self.lambda)/(2*sin(self.phi))^(2*self.lambda);
