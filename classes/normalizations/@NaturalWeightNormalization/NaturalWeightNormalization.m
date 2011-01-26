classdef NaturalWeightNormalization < WeightNormalization
% NaturalWeightNormalzation -- Specifies a "natural" weight normalization
%
% self = NaturalWeightNormalization()
%
%     A singleton class serving as a fingerprint indicating that a weight
%     function has a "natural" normalization.  "Natural" means the
%     normalization that is prevalent in the literature, with the caveat that if
%     an affine map results in a non-standard domain, the Jacobian of the affine
%     map is built into the weight function so that the e.g. "orthonormal"
%     evaluations are not changed.

  properties(Access=private)
    description;
  end

  methods(Access=private)
    function self = NaturalWeightNormalization()
      self.description = 'Natural weight normalization';
    end
  end

  methods(Static)
    function obj = instance()
      persistent self
      if isempty(self)
        obj = NaturalWeightNormalization();
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
