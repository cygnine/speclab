function[composition] = compose(self,other)
% compose -- Composition of affine maps
%
% composition = compose(self, other)
%
%     Given another AffineMap instance "other", this function creates the
%     composition map f = self \circ other.

% Check for dimensions:
if self.n ~= other.m
  error('Affine map domains and ranges do not have commensurate dimensions');
end

composition = AffineMap(self.A*other.A, self.b + self.A*other.b);
