function[C] = coefficient_tensorization(self,coefficients)
% coefficient_tensorization -- tensorizes univariate expansion coefficients
%
% C = coefficient_tensorization(self, coefficients)
%
%     The input coefficients is a length-(self.dim) cell array, with entry
%     d being an N x M(d) array. Row n of coefficients{d} represents expansion
%     coefficients of a function in the univariate basis self.bases{d}. The
%     columns are ordered by the univariate indexing for self.bases{d}.
%
%     C is an N x R array, where row n of C is expansion coefficients of the
%     function that is the product of coefficients{d} for all d. The number of
%     columns R is chosen to be as small as possible. C is returned as a sparse
%     array.

if length(coefficients) ~= self.dim
  error('Input cell array must have self.dim elements');
end

N = size(coefficients{1},1);
for d = 2:self.dim
  if size(coefficients{d},1) ~= N
    error('Input matrices in cell array must have the same number of rows');
  end
end

% Attempt to precompute total degree we'll need
maxN = zeros([self.dim 1]);
for d = 1:self.dim
  maxN = size(coefficients{d},2) - 1; % poly degree in dimension d
end

R = 
