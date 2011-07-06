function[result] = mtimes(self, other)
% mtimes -- Multiply together two FactoredPolynomial's
%
% result = mtimes(self, other)
%
%     If 'other' is a FactoredPolynomial instance, then result is the
%     FactoredPolynomial object that is the multiplication of self and other.
%
%     If 'other' is just a scalar numeric, then result is just self with a
%     modified leading_coefficient.

if isa(other, 'numeric')
  result = self;
  self.leading_coefficient = self.leading_coefficient*other;
elseif isa(other, 'FactoredPolynomial')
  result = self;
  result.roots = [self.roots; other.roots];
  result.leading_coefficient = self.leading_coefficient*other.leading_coefficient;
else
  error(['Cannot multiply together object of type ' class(self) ' with type ' class(other)]);
end
