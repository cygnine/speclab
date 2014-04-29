function w = weight(self, x, d)
% weight -- Weight function for Hermite polynomials
%
% w = weight(self, x, d)
% 
%     On the standard interval [-Inf,Inf], define the weight function
%
%       w(r) = |r|^(2*self.mu) * exp(-r^2)
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
  w = abs(r).^(2*self.mu).*exp(-r.^2);
case 1
  if self.mu==0
    w = -2*r.*exp(-r.^2);
  else
    w = 2*exp(-r.^2).*(-r.*abs(r).^(2*self.mu) + 2*self.mu*sgn(r).*abs(r).^(2*self.mu-1));
  end
  % Jacobian factor for derivative:
  w = w*self.map_to_standard_domain.A;
end

w = self.scale_weight(w);
