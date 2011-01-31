function w = scale_weight(self,w)
% w = scale_weight(self,w)
%
%     Assumes that the weight function input w has 'default' weight scaling
%     (the definition of which is dependent on the particular polynomial
%     family). This function then rescales w to match the custom weight given by
%     self.weight_normalization.

persistent classical natural probability
if isempty(classical)
  classical = ClassicalWeightNormalization.instance();
  natural = NaturalWeightNormalization.instance();
  probability = ProbabilityWeightNormalization.instance();
end

if self.weight_normalization==classical
  % Nothing happens
  return
elseif self.weight_normalization==natural
  % We place the Jacobian as part of the weight function
  w = self.map_to_standard_domain.A*w;
elseif self.weight_normalization==probability
  % Just divide by the integral on the standard interval
  [a,b] = self.recurrence(0);
  w = self.map_to_standard_domain.A*w/b;
else
  error('This basis does not support the given function normalization type');
end
