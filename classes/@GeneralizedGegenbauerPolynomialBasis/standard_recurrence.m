function[a,b] = standard_recurrence(self, n)
% [a,b] = standard_recurrence(self, n)
%
%     Returns the recurrence coefficients with non-negative indices n for the
%     Generalized Gegenbauer polynomials defined by (self.lambda, self.mu). 

N = max(n) + 1;

[lambda,mu] = deal(self.lambda, self.mu);
%[alpha,beta] = deal(self.alpha(dim), self.beta(dim));
a = zeros(size(n));
b = ones(size(n));

% Initial conditions:
flags0 = (n==0);
if any(flags0)
  b(flags0) = exp(- gammaln(lambda+mu+1) + gammaln(mu+1/2) + gammaln(lambda+1/2));
end

flags1 = (n==1);
b(flags1) = (2*mu+1)/(2*(lambda+mu+1));

flags = not(flags0 | flags1);
flags_odd = flags & (mod(n,2)==1);
flags_even = flags & (mod(n,2)==0);
b(flags_odd) = (n(flags_odd) + 2*mu).*(n(flags_odd) + 2*lambda - 1 + 2*mu)./...
                (4.*(n(flags_odd) + lambda + mu - 1).*(n(flags_odd) + lambda + mu));
b(flags_even) = n(flags_even).*(n(flags_even) + 2*lambda - 1)./...
                (4.*(n(flags_even) + lambda + mu - 1).*(n(flags_even) + lambda + mu));
