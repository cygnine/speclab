function w = weight(self, x, d)
% weight -- Weight function for Jacobi polynomials
%
% w = weight(self, x, d)
% 
%     On the standard interval [-1,1], define the weight function
%
%       w(r) = [ (1-r)^self.alpha ] x [ (1+r)^self.beta ].
%
%     On the physical domain x(r), this function evaluates the d'th derivative
%     of C*w(r(x)), where C is determined by the scaling in self.scale_weight
%     and self.weight_normalization.

if nargin < 3
  d = 0;
end

r = self.map_to_standard_domain(x);

switch d
case 0
  w = ((1-r).^self.alpha).*((1+r).^self.beta);
case 1
  w = ((1-r).^(self.alpha-1)).*((1+r).^(self.beta-1));
  if self.alpha==0
    w1 = 0;
  else
    w1 = w.*(1+r)*self.alpha;
  end
  if self.beta==0
    w2 = 0;
  else
    w2 = w.*(1-r)*self.beta;
  end
  w = w1 + w2;
end

w = self.scale_weight(w);
