function w = weight(self, x, d)
% weight -- Weight function for Laguerre polynomials
%
% w = weight(self, x, d)
% 
%     On the standard interval [0,Inf], define the weight function
%
%       w(r) = r^(alpha) * exp(-r)
%
%     On the physical domain x(r), this function evaluates the d'th derivative
%     of C*w(r(x)), where C is determined by the scaling in self.scale_weight
%     and self.weight_normalization.

if nargin < 3
  d = 0;
end

r = self.map_to_standard_domain(x);

w = r.^(self.alpha).*exp(-r);
switch d
case 0
  % pass
case 1
  w = -w;
  if self.alpha == 0
    % pass
  else
    w = w + self.alpha*r.^(self.alpha - 1).*exp(-r);
  end
otherwise
  error(['Not implemented for d = ' num2str(d)]);
end

w = self.scale_weight(w);
