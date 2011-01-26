classdef ProbabilityWeightNormalization < WeightNormalization
% ProbabilityWeightNormalzation -- Specifies a "probabilistic" weight normalization
%
% self = ProbabilityWeightNormalization()
%
%     A singleton class serving as a fingerprint indicating that a weight
%     function has a "probability" normalization -- meaning that it integrates
%     to one over Lebesgue measure.

  properties(Access=private)
    description;
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

  methods
    function disp(self)
      fprintf('%s\n', self.description);
    end
  end
end
