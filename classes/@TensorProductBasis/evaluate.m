function[V] = evaluate(self,x,n)
% evaluate -- Evaluates tensor-product basis
%
% V = evaluate(self, x, n)
%
%     Evaluates the (one-based-indexed) n'th basis functions at the locations x.
%     x is an array of size N x self.dim, where each row is a point in
%     d-dimensional space. The output V is an array of size N x length(n), where
%     each column represents the appropriate basis function n evaluated at all
%     locations x.

n = n(:);
if size(x,2) ~= self.dim
  error('Input points must have the same dimension (columns) as self.dim');
end
if any(n<1)
  error('Tensor-product bases have linear indexing that is one-based');
end

inds = self.indexing(n-1, 'dim', self.dim);

V = ones(size(x,1), length(n));

for q = 1:self.dim
  V = V.*self.bases{q}(x(:,q), inds(:,q));
end
