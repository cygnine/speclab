function[a,b] = recurrence(self, n)
% [a,b] = recurrence(n,d=1);
%
%     Returns the recurrence coefficients with indices n for the Jacobi
%     polynomials defined by (alpha, beta). 

%if nargin < 3
%  dim = 1;
%end

N = max(n) + 1;
[alpha,beta] = deal(self.alpha, self.beta);
%[alpha,beta] = deal(self.alpha(dim), self.beta(dim));
a = (beta^2-alpha^2)*ones(size(n));
b = ones(size(n));

% Initial conditions:
flags0 = (n==0);
if any(flags0)
  a(flags0) = (beta-alpha)/(alpha+beta+2);
  b(flags0) = 2^(alpha+beta+1)*gamma(alpha+1)*gamma(beta+1)/gamma(alpha+beta+2);
end

flags1 = (n==1);
a(flags1) = a(flags1)./((2+alpha+beta)*(4+alpha+beta));
b(flags1) = 4*(1+alpha)*(1+beta)/((2+alpha+beta)^2*(3+alpha+beta)); 

flags = not(flags0 | flags1);
a(flags) = a(flags)./((2*n(flags)+alpha+beta).*(2*n(flags)+alpha+beta+2));
b(flags) = 4*n(flags).*(n(flags)+alpha).*(n(flags)+beta).*(n(flags)+alpha+beta);
b(flags) = b(flags)./((2*n(flags)+alpha+beta).^2.*...
                     (2*n(flags)+alpha+beta+1).*(2*n(flags)+alpha+beta-1));
