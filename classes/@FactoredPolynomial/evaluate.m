function[p] = evaluate(self, x)
% evaluate -- Evaluates the FactoredPolynomial
%
% p = evaluate(self,x)
%
%     Evaluates the FactoredPolynomial at the locations x.

p = self.leading_coefficient*ones(size(x));
for q = 1:self.degree
  p = p.*(x - self.roots(q));
end
