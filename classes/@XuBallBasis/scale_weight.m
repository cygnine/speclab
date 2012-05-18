function w = scale_weight(self,w)
% w = scale_weight(self,w)
%
%     Assumes that the weight function input w has 'default' weight scaling
%     (the definition of which is dependent on the particular polynomial
%     family). This function then rescales w to match the custom weight given by
%     self.weight_normalization.

if self.weight_normalization=='probability'
  % Nothing happens
  w = w*(2*self.mu+1)/(2*pi)*det(self.map_to_standard_domain.A);
else
  w = scale_weight@MultivariateOrthogonalPolynomialBasis(self, w);
end
