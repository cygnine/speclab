classdef ProbabilistFunctionNormalization < FunctionNormalization
% ProbabilistFunctionNormalzation -- Specifies a "probabilist" function normalization
%
% self = ProbabilistFunctionNormalization()
%
%     A singleton class serving as a fingerprint indicating that a basis set has
%     function evaluations that correspond to a "probabilist" normalization.
%     "Claisscal" means the normalization that is prevalent in the literature.

  properties(SetAccess=private)
    description;
    ids = {'probability', 'prob', 'probabilist'};
  end

  methods(Access=private)
    function self = ProbabilistFunctionNormalization()
      self.description = 'Probabilist function normalization';
    end
  end

  methods(Static)
    function obj = instance()
      persistent self
      if isempty(self)
        obj = ProbabilistFunctionNormalization();
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
