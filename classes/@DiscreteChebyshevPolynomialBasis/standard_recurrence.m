function[a,b] = standard_recurrence(self, n)
% [a,b] = standard_recurrence(self, n)
%
%     Returns the recurrence coefficients with non-negative indices n for the
%     discrete Chebyshev polynomials defined by self.N.

N = max(n);
assert(N <= self.N-1, ['Indexing for this basis is defined only up to ' num2str(self.N)]);

nzero = (n==0);

a = self.N/2*(1 - 1/self.N)*ones(size(n));

b = zeros(size(a));
b(nzero) = self.N;
% And now the non-zero terms
nzero = ~nzero;
b(nzero) = self.N^2/4*(1 - (n(nzero)/self.N).^2)./(4 - 1./n(nzero).^2);
