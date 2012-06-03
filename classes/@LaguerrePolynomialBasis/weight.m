function w = weight(self, x)
% weight -- Weight function for Laguerre polynomials
%
% w = weight(self, x)
% 
%     On the standard interval [0,Inf], define the weight function
%
%       w(r) = r^(alpha) * exp(-r)
%
%     On the physical domain x(r), this function evaluates C*w(r(x)), where C is
%     determined by the scaling in self.scale_weight and
%     self.weight_normalization.

r = self.map_to_standard_domain(x);
w = r.^(self.alpha).*exp(-r);
w = self.scale_weight(w);
