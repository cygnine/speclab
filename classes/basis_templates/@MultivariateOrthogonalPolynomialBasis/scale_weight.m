function w = scale_weight(self,w)
% w = scale_weight(self,w)
%
%     Assumes that the weight function input w has 'default' weight scaling
%     (the definition of which is dependent on the particular polynomial
%     family). This function then rescales w to match the custom weight given by
%     self.weight_normalization.

if self.weight_normalization=='classical'
  % Nothing happens
  return
elseif self.weight_normalization=='natural'
  % We place the Jacobian as part of the weight function
  w = det(self.map_to_standard_domain.A)*w;
else
  error('This basis does not support the given function normalization type');
end
