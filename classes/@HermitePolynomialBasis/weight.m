function w = weight(self, x)
% weight -- Weight function for Hermite polynomials
%
% w = weight(self, x)
% 
%     On the standard interval [-Inf,Inf], define the weight function
%
%       w(r) = |r|^(2*self.mu) * exp(-r^2)
%
%     On the physical domain x(r), this function evaluates C*w(r(x)), where C is
%     determined by the scaling in self.scale_weight and
%     self.weight_normalization.

r = self.map_to_standard_domain(x);
w = abs(r).^(2*mu)*exp(-r.^2);
w = self.scale_weight(w);
