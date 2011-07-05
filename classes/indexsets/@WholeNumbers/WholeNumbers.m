classdef WholeNumbers < IndexSet

  properties(SetAccess=private)
    descriptor
  end
  methods(Static)
    function self = instance()
      persistent obj
      if isempty(obj)
        obj = WholeNumbers();
      end
      self = obj;
    end
  end
  methods(Access=private)
    function self = WholeNumbers()
      % WholeNumbers -- The whole numbers 0, 1, 2, ...
      %
      % self = WholeNumbers()
      %
      %     The set of whole numbers (non-negative integers). This is a class
      %     with a singleton instance.
      self.descriptor = 'Whole numbers';
    end
  end
  methods
    validate(self,inp)
  end
end
