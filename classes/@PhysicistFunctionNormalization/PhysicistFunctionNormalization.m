classdef PhysicistFunctionNormalization < FunctionNormalization
% PhysicistFunctionNormalzation -- Specifies a "physicist" function normalization
%
% self = PhysicistFunctionNormalization()
%
%     A singleton class serving as a fingerprint indicating that a basis set has
%     function evaluations that correspond to a "physicist" normalization.
%     "Claisscal" means the normalization that is prevalent in the literature.

  properties(Access=private)
    description;
  end

  methods(Access=private)
    function self = PhysicistFunctionNormalization()
      self.description = 'Physicist function normalization';
    end
  end

  methods(Static)
    function obj = instance()
      persistent self
      if isempty(self)
        obj = PhysicistFunctionNormalization();
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
