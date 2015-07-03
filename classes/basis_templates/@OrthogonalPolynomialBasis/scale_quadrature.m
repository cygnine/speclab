function[x,w] = scale_quadrature(self,x,w)
% [x,w] = scale_quadrature(self,x,w)
%
%     This function scales quadrature rules living on standard intervals with
%     standard weight normalization to 

persistent classical natural probability
if isempty(classical)
  classical = ClassicalWeightNormalization.instance();
  natural = NaturalWeightNormalization.instance();
  probability = ProbabilityWeightNormalization.instance();
end

% No matter what:
x = self.map_to_domain(x);

if self.weight_normalization==classical
  w = w*self.map_to_domain.A;
  % Nothing happens
  return
elseif self.weight_normalization==natural
  % We place the Jacobian as part of the weight function
  % w = self.map_to_standard_domain.A*w;
  % do nothing
elseif self.weight_normalization==probability
  % Just divide by the integral on the standard interval
  [a,b] = self.standard_recurrence(0);
  %w = self.map_to_standard_domain.A*w/b;
  w = w/b;
else
  error('This basis does not support the given function normalization type');
end
