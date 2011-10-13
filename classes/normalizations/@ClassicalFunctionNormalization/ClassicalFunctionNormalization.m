classdef ClassicalFunctionNormalization < FunctionNormalization
% ClassicalFunctionNormalzation -- Specifies a "classical" function normalization
%
% self = ClassicalFunctionNormalization()
%
%     A singleton class serving as a fingerprint indicating that a basis set has
%     function evaluations that correspond to a "classical" normalization.
%     "Claisscal" means the normalization that is prevalent in the literature.
%
%     Use the syntax classname.instance() to instantiate this class.

  properties(SetAccess=protected)
    ids = {'classical', 'class'};
  end

  methods(Access=private)
    function self = ClassicalFunctionNormalization()
      self.description = 'Classical function normalization';
    end
  end

  methods(Static)
    function obj = instance()
      persistent self
      if isempty(self)
        obj = ClassicalFunctionNormalization();
        self = obj;
      else
        obj = self;
      end
    end
  end

end
