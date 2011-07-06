function[p] = evaluate(self, x)
% evaluate -- Evaluates the MeasureModificationPolynomial
%
% p = evaluate(self,x)
%
%     Evaluates the MeasureModificationPolynomial at the locations x.

p = self.leading_coefficient*ones(size(x));
for q = 1:self.num_exterior
  p = p.*(x - self.exterior_roots(q));
end
for q = 1:self.num_conjugate
  p = p.*(x.^2 + abs(self.conjugate_roots(q)).^2);
end
for q = 1:self.num_quadratic
  p = p.*(x - self.quadratic_roots(q)).^2;
end
