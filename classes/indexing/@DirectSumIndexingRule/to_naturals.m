function[n] = to_naturals(self, a)
% to_naturals -- mapping to index set
%
% n = to_naturals(self, a)
%
%     a is an N x self.dim array of indices from this rules, and n is then a
%     length-N column vector containing the corresponding natural indices.

N = size(a,1);
if size(a,2) ~= self.dim
  error('Number of columns of input does not match dimension of this rule');
end

for q = 1:self.dim
  a(:,q) = self.rules{q}.inv(a(:,q));
end
n = self.multirule.inv(a);
