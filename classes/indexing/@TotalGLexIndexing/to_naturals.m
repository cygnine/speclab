function[n] = to_naturals(self, a)
% to_naturals -- transforms array indexing to linear indexing
%
% [n] = to_naturals(self, a)
%
%      Computes one-based indexing for lexicographic ordering on
%      the self.dim-dimensional Cartesian product of N_0, where N_0 is the
%      naturals unioned with 0.

persistent subdim spdim
if isempty(subdim)
  from speclab.common.tensor import subspace_dimension as subdim
  from speclab.common.tensor import space_dimension as spdim
end

if size(a,1) ~= self.dim
  error('The input must have row size equal to self.dim');
end

a = a.';
tempsum = fliplr(cumsum(fliplr(a),2));

% Copy-paste from TotalReverseLexicographicIndexing
n = zeros([size(a,1) 1]);
for q = 1:self.dim
  n = n + spdim(self.dim-q+1, tempsum(:,q)-1);
end
n = n + 1;

nds = spdim(self.dim, tempsum(:,1)-1);

% "Reverse" the intra-degree ordering
%n = nds + (subdim(self.dim, tempsum(:,1)) + 1 - (n - nds));
% The above is the same as:
n = 1 + 2*nds - n + subdim(self.dim, tempsum(:,1));
