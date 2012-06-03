function[V] = evaluate(self, x, n)
% p = evaluate(x,n)
%
%     Evaluates the orthogonal polynomial basis at the locations x on the 2D
%     disc. x must be a vector with two columns (2D Cartesian coordinates). p
%     is a size(x,1) x N array, where N is the number of polynomials requested
%     by the index set n.

[n_array, nsize, numeln] = self.indexing(n);

% Find maximum degree, and generate needed x-bases.
total_degree = sum(n_array,2);
N = max(sum(total_degree,2));
k = length(self.univariate_bases{1});
while k < N+1
  self.univariate_bases{1}{end+1} = JacobiPolynomialBasis('alpha', self.mu + k, 'beta', self.mu + k);
  k = k + 1;
end

Vsize = [size(x,1), numeln];
V = ones(Vsize);

% First generate x contributions: we have to do this degree-by-degree
K = max(n_array(:,2));
for k = 0:K
  kflags = find(n_array(:,2)==k);
  if not(isempty(kflags))
    V(:,kflags) = V(:,kflags).*self.univariate_bases{1}{k+1}(x(:,1), n_array(kflags,1));
  end
end

% Now generate y contributions: we have to be careful at x= \pm 1 for this
xflags = (abs(x(:,1)) ~= 1);
xweight = repmat(1-x(xflags,1).^2, [1 numeln]);
xweight = xweight.^repmat(n_array(:,2).'/2, [sum(xflags) 1]);
V(xflags,:) = V(xflags,:).*xweight;
V(xflags,:) = V(xflags,:).*self.univariate_bases{2}(x(xflags,2)./sqrt(1 - x(xflags,1).^2), n_array(:,2));

xflags = ~xflags;
% At boundaries, just multiply by leading coefficient contribution of the y-Gegenbauer contribution 
V(xflags,:) = V(xflags,:).*...
        [repmat(x(xflags,2), [1 numeln]).^repmat(n_array(:,2).', [sum(xflags) 1])]*...
         diag(self.univariate_bases{2}.leading_coefficient(n_array(:,2)));

V = self.scale_functions(V,n_array);

if (size(V,3)==1) && (size(x,1)==1)
  V = reshape(V, nsize);
end
