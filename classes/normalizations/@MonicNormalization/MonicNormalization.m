classdef MonicNormalization < FunctionNormalization
% MonicNormalization -- Specifies a "monic" function normalization
%
% self = MonicNormalization()
%
%     A singleton class serving as a fingerprint indicating that a basis set has
%     function evaluations that correspond to a "monic" normalization -- i.e.
%     the leading term has coefficient 1.

  properties(SetAccess=private)
    description;
    ids = {'monic', 'mon'};
  end

  methods(Access=private)
    function self = MonicNormalization()
      self.description = 'Monic function normalization';
    end
  end

  methods(Static)
    function obj = instance()
      persistent self
      if isempty(self)
        obj = MonicNormalization();
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
