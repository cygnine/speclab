classdef WholeNumbers < IndexSet
% WholeNumbers -- The set of natural numbers
%
% self = WholeNumbers()
%
%     An IndexSet singleon subclass corresponding to the natural numbers.
%
% WholeNumbers Properties:
%   descriptor - A human-readable description: 'Whole numbers'
% WholeNumbers Methods:
%   validate - Determines is an input is from the index set
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
