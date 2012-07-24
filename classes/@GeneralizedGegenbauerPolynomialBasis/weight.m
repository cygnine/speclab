function w = weight(self, x)
% weight -- Weight function for Generalized Gegenbauer polynomials
%
% w = weight(self, x)
% 
%     On the standard interval [-1,1], define the weight function
%
%       w(r) =  (1-r^2)^self.alpha x |r|^(2*self.mu)
%
%     On the physical domain x(r), this function evaluates C*w(r(x)), where C is
%     determined by the scaling in self.scale_weight and
%     self.weight_normalization.

r = self.map_to_standard_domain(x);
w = ((1-r.^2).^(self.lambda-1/2)).*(abs(r).^(2*self.mu));
w = self.scale_weight(w);
