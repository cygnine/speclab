function[lmbdaout] = derivative_expansion(self, n, d)
% derivative_expansion -- expansion coefficients of polynomial derivatives
%
% lmbda = derivative_expansion(self, n, d)
%
%     n is a vector of valid indices, and d is a vector of positive integers.
%
%     If length(d)==1: returns a length(n) x M matrix with expansion
%     coefficients of self such that the d'th derivative of the n'th orthogonal
%     polynomials has expansion coefficients lmbda. Row q+1 of lmbda has
%     coefficients:
%
%         p^{(d)}_{n(q)} = \sum_{m=0}^{M-1} lmbda(q,m+1) * p_{m},
%
%     where for clarity above we have assumed that n and m are
%     ZeroBasedIndexing. p^{(d)}_n is the d'th derivative of the orthogonal
%     polynomial p_n associated with self. The number of columns M of lmbda is
%     chosen to be as small as possible: M = max(n+1) - d.
%
%     If d is a vector: returns a three-dimensional array lmbda, where
%     lmbda(:,:,r) is the matrix lmbda corresponding to d = d(r). The number of
%     columns of lmbda is the numbers of columns for the smallest value of d.

persistent spdiag
if isempty(spdiag)
  from labtools import spdiag
end

% First change indices to ZeroBasedIndexing:
[n, nsize, numeln] = self.indexing(n);
N = max(n);

d = d(:);
maxd = max(d);

M = N-1;
% First let's code for length(d)==1
%lmbda = zeros([length(n) N+1-d]);
lmbdaout = zeros([N+1 M+1 length(d)]);
prevlmbda = speye([N+1 M+1]);

% Need trivial output?
inds  = find(d==0);
for q = 1:length(inds)
  lmbdaout(:,:,q) = prevlmbda;
end

[a,b] = self.recurrence(0:N);
b = sqrt(b);

% Iterate over derivatives
for D = 1:maxd
  % First initial conditions:
  lmbda = zeros([N+1 M+1]);

  % lmbda^{(D)}_{D,0}
  lmbda(D+1,1) = factorial(D)/prod(b(2:D+1));

  % And now iteration
  % lmbda^{(D)}_{ns,:}
  for ns = D+1:N
    ms = 0:(ns-D);
    lmbda(ns+1,ms+1) = lmbda(ns+1,ms+1) - b(ns)*lmbda(ns-1,ms+1) + ...
                       (a(ms+1) - a(ns)).*lmbda(ns,ms+1) + ...
                       D*prevlmbda(ns,ms+1);

    %lmbda(ns+1,m+1) = lmbda(ns+1,m+1) + b(m+2)*lmbda(ns,m+2);
    lmbda(ns+1,1:M) = lmbda(ns+1,1:M) + b(2:M+1).*lmbda(ns,2:M+1);

    %lmbda(ns+1,m+1) = lmbda(ns+1,m+1) + b(m+1)*lmbda(ns,m);
    lmbda(ns+1,2:M+1) = lmbda(ns+1,2:M+1) + b(2:M+1).*lmbda(ns,1:M);

    % lmbda^{(D)}_{ns,m}
    lmbda(ns+1,:) = lmbda(ns+1,:)/b(ns+1);
  end
  prevlmbda = lmbda;

  % Need any output?
  inds = find(d==D);
  for q = 1:length(inds)
    lmbdaout(:,:,inds(q)) = lmbda;
  end
end

d = min(d);
lmbdaout = lmbdaout(n+1,1:(N+1-d),:);

% And now scale functions:
colscale = spdiag(1./self.scale_functions(ones([1 M+1]), 0:M));
rowscale = spdiag(self.scale_functions(ones([1 length(n)]), n));
for D = 1:size(lmbdaout,3)
  lmbdaout(:,:,D) = rowscale*lmbdaout(:,:,D)*colscale;
end
