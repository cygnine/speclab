classdef IndexSet < Singleton
% IndexSet -- A superclass for index sets
%
% self = IndexSet()
%
%     A superclass for all singleton classes that represent index sets.  There
%     are no public data properties for this class, and this class has no
%     constructor.
  properties(SetAccess=private,Abstract=true)
    descriptor
  end
  methods
    function disp(self)
      fprintf('%s set of indices\n', self.descriptor);
    end
  end
  methods(Abstract=true)
    validate
  end
end
