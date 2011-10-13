classdef ProbabilityWeightNormalization < WeightNormalization
% ProbabilityWeightNormalzation -- Specifies a "probabilistic" weight normalization
%
% self = ProbabilityWeightNormalization()
%
%     A singleton class serving as a fingerprint indicating that a weight
%     function has a "probability" normalization -- meaning that it integrates
%     to one over Lebesgue measure.
%
%     Use the syntax classname.instance() to instantiate this class.

  properties(SetAccess=protected)
    ids = {'probability', 'prob'};
  end

  methods(Access=private)
    function self = ProbabilityWeightNormalization()
      self.description = 'Probability weight normalization';
    end
  end

  methods(Static)
    function obj = instance()
      persistent self
      if isempty(self)
        obj = ProbabilityWeightNormalization();
        self = obj;
      else
        obj = self;
      end
    end
  end

end
