function w = weight(self, x)
% weight -- Weight function for Jacobi polynomials
%
% w = weight(self, x)
% 
%     On the standard interval [-1,1], define the weight function
%
%       w(r) = [ (1-r)^self.alpha ] x [ (1+r)^self.beta ].
%
%     On the physical domain x(r), this function evaluates C*w(r(x)), where C is
%     determined by the scaling in self.scale_weight and
%     self.weight_normalization.

r = self.map_to_standard_domain(x);
w = ((1-r).^self.alpha).*((1+r).^self.beta);
w = self.scale_weight(w);
