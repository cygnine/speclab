classdef IntegerNumbers < IndexSet

  properties(SetAccess=private)
    descriptor
  end
  methods(Static)
    function self = instance()
      persistent obj
      if isempty(obj)
        obj = IntegerNumbers();
      end
      self = obj;
    end
  end
  methods(Access=private)
    function self = IntegerNumbers()
      % IntegerNumbers -- The integers ... -3, -2, -1, 0, 1, 2, 3, ...
      %
      % self = IntegerNumbers()
      %
      %     The set of integers. This is a class with a singleton instance.
      self.descriptor = 'Integer numbers';
    end
  end
  methods
    validate(self,inp)
  end
end
