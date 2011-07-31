classdef NaturalNumbers < IndexSet
% NaturalNumbers -- The set of natural numbers
%
% self = NaturalNumbers()
%
%     An IndexSet singleon subclass corresponding to the natural numbers.
%
% NaturalNumbers Properties:
%   descriptor - A human-readable description: 'Natural numbers'
% NaturalNumbers Methods:
%   validate - Determines is an input is from the index set
  properties(SetAccess=private)
    descriptor
  end
  methods(Static)
    function self = instance()
      persistent obj
      if isempty(obj)
        obj = NaturalNumbers();
      end
      self = obj;
    end
  end
  methods(Access=private)
    function self = NaturalNumbers()
      % NaturalNumbers -- The natural numbers 1, 2, 3, ...
      %
      % self = NaturalNumbers()
      %
      %     The set of natural numbers (positive integers). This is a class with
      %     a singleton instance.
      self.descriptor = 'Natural numbers';
    end
  end
  methods
    validate(self,inp)
  end
end
