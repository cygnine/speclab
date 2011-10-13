classdef ClassicalWeightNormalization < WeightNormalization
% ClassicalWeightNormalzation -- Specifies a "classical" weight normalization
%
% self = ClassicalWeightNormalization()
%
%     A singleton class serving as a fingerprint indicating that a weight
%     function has a "classical" normalization.  "Claisscal" means the
%     normalization that is prevalent in the literature -- this can be construed
%     as a default normalization. When an affine map is present, the weight
%     function evaluates to the value of the mapped coordinate (any resulting
%     Jacobian is built into the function evaluations).
%
%     Use the syntax classname.instance() to instantiate this method.

  properties(SetAccess=protected)
    ids = {'classical', 'class'};
  end

  methods(Access=private)
    function self = ClassicalWeightNormalization()
      self.description = 'Classical weight normalization';
    end
  end

  methods(Static)
    function obj = instance()
      persistent self
      if isempty(self)
        obj = ClassicalWeightNormalization();
        self = obj;
      else
        obj = self;
      end
    end
  end

end
