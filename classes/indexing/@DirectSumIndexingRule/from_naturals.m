function[a] = from_naturals(self, n)
% from_naturals -- mapping to index set
%
% a = from_naturals(self, n)
%
%     n is an array of natural numbers, and the resultant array a has rows equal
%     to numel(n), and columns equal to self.dim.

a = self.multirule(n(:));
for q = 1:self.dim
  a(:,q) = self.rules{q}(a(:,q));
end
