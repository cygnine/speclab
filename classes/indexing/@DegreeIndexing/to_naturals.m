function[n] = to_naturals(self, a)
% to_naturals -- transforms degree indexing to linear indexing
%
% [n] = to_naturals(self, a)
%
%      Computes one-based indexing for degree-based ordering on the whole numbers.

if all(size(a) ~= 1)
  warning('Input should be a vector. Flattening array to vector');
  a = a(:);
end
if any(a < 0)
  error('Input must be a non-negative integer');
end

N = max(a);

n = [];
sizes_minus1 = self.total_polynomial_space_dimension(self.dim, a-1);
sizes = self.total_polynomial_space_dimension(self.dim, a);
for q = 1:length(a)
  current_contribution = (sizes_minus1(q)+1):sizes(q);
  n = [n; current_contribution(:)];
end

% Take care of case where a is a row vector
if size(a,2) > 1
  n = n.';
end
