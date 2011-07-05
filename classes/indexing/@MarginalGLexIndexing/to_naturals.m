
% "Reverse" the intra-degree ordering
%n = nds + (subdim(self.dim, tempsum(:,1)) + 1 - (n - nds));
% The above is the same as:
n = 1 + 2*nds - n + subdim(self.dim, tempsum(:,1));

function[n] = to_naturals(self, a)
% to_naturals -- transforms array indexing to linear indexing
%
% [n] = to_naturals(self, a)
%
%      Computes one-based indexing for lexicographic ordering on
%      the self.dim-dimensional Cartesian product of N_0, where N_0 is the
%      naturals unioned with 0.

persistent spdim subdim
if isempty(spdim)
  %from speclab.common.tensor import space_dimension as spdim
  from speclab.common.tensor import marginalpoly_space_dim as spdim
  from speclab.common.tensor import marginalpoly_subspace_dim as subdim
end

if size(a,1) ~= self.dim
  error('The input must have row size equal to self.dim');
end
a = a.';

% Copy+paste from MarginalGRevLexIndexing.to_naturals

% First add in marginal degree:
[marginal_degrees, columns] = max(a, [], 2);
n = spdim(self.dim, marginal_degrees-1);

% Loop over each column where marginal appears:
for q = 1:self.dim
  flags = (columns==q);

  % Add in 'prefix' terms: shift by dimensions of 'previous' subspaces
  if q > 1
    for qq = 1:(q-1);
      factors = marginal_degrees(flags) - a(flags,qq);
      n(flags) = n(flags) + (factors>0).*spdim(self.dim-qq, marginal_degrees(flags));
      factors = max(zeros(size(factors)), factors-1);
      n(flags) = n(flags) + factors.*subdim(self.dim-qq, marginal_degrees(flags));
    end
  end

  % And now add in suffix terms: factor in `floating point' remainder
  if q < self.dim
    % *sigh* sub2ind not vectorized in 'SIZ' input
    base = marginal_degrees(flags)+1;
    for qq = 1:(self.dim-q)
      significand = marginal_degrees(flags) - a(flags,q+qq);
      exponent = self.dim - q - qq;
      n(flags) = n(flags) + significand.*base.^exponent;
    end
  end
end

% Now reverse:
%n = n - spdim(self.dim, marginal_degrees-1);
%n = subdim(self.dim, marginal_degrees) - n;
%n = n + spdim(self.dim, marginal_degrees-1);

% Save as above:
n = 2*spdim(self.dim, marginal_degrees-1) + subdim(self.dim, marginal_degrees) - n;
