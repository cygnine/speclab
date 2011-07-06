classdef FactoredPolynomial
% FactoredPolynomial -- Univariate polynomial in factored form
%
% self = FactoredPolynomial(roots,  [[ leading_coefficient = 1 ]] )
%
%     Object representation of a scalar-valued univariate polynomial p in
%     factored form:
%
%      p(x) = leading_coefficient * (x - roots(1)) * (x - roots(2)) * ...
%
%     Factored representations are useful in certain situations; otherwise
%     one might as well use Matlab's builtin representation using a row
%     vector of monomial coefficients.
  properties
    roots
    degree
    leading_coefficient = 1
    coefficients
  end
  methods
    function[self] = FactoredPolynomial(roots, leading_coefficient)
      if nargin < 2
        leading_coefficient = 1;
      end

      if not(isnumeric(roots)) || not(isnumeric(leading_coefficient))
        error('Input roots and leading coefficient must be numeric values')
      end
      if numel(leading_coefficient) > 1
        warning('Leading coefficient must be a scalar -- extracting first element of array and discarding the rest'); 
      end

      self.roots = roots(:);
      self.leading_coefficient = leading_coefficient(1);
    
    end
    function[self] = set.roots(self,newroots)
      % Update degree, set new roots, and recompute monomial coefficients
      self.degree = length(newroots);
      self.roots = newroots;
      self.coefficients = 1;
      for q = 1:self.degree
        self.coefficients = conv(self.coefficients, [1 -self.roots(q)]);
      end
    end
    function[self] = set.leading_coefficient(self,newcoeff)
      % Update monomial coefficients and set new leading coefficient
      self.coefficients = self.coefficients/self.leading_coefficient;
      self.coefficients = self.coefficients*newcoeff;
      self.leading_coefficient = newcoeff;
    end
    p = evaluate(self,x)
    result = mtimes(self,other)
  end
end
