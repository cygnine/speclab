classdef OrthonormalNormalization < FunctionNormalization
% OrthonormalNormalization -- Specified an "orthonormal" function normalization
%
% self = OrthonormalNormalization()
%
%     A singleton class serving as a fingerprint indicating that a basis set has
%     function evaluations that correspond to an "orthonormal" normalization.
%     Here "orthonormal" refers to unit norm under some innate norm (usually an
%     L^2 inner product) that is defined by the basis set.

  properties(SetAccess=private)
    description;
    ids = {'orthonormal', 'normal'};
  end

  methods(Access=private)
    function self = OrthonormalNormalization()
      self.description = 'Orthonormal function normalization';
    end
  end

  methods(Static)
    function obj = instance()
      persistent self
      if isempty(self)
        obj = OrthonormalNormalization();
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
