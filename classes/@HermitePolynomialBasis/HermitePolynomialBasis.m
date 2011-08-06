classdef HermitePolynomialBasis < OrthogonalPolynomialBasis
  properties
    mu
  end
  methods
    function self = HermitePolynomialBasis(varargin)
    % self = HermitePolynomialBasis(varargin)
    %
    %     Creates an instance of a Hermite Polynomial spectral basis.

      persistent inparse
      if isempty(inparse)
        inparse = inputParser();
        inparse.KeepUnmatched=true;

        inparse.addParamValue('mu', 0);
        inparse.addParamValue('normalization', 'normal');
      end
      %inputs = {'mu', 'normalization'};
      %defaults = {0, 'normal'};

      inparse.parse(varargin{:});
      opt = inparse.Results;

      self = self@OrthogonalPolynomialBasis(varargin{:});

      self.allowed_function_normalizations{end+1} = ProbabilistFunctionNormalization.instance();
      self.allowed_function_normalizations{end+1} = PhysicistFunctionNormalization.instance();
      self.allowed_function_normalizations{end+1} = ClassicalFunctionNormalization.instance();
      self.normalization = self.function_normalization_parser(opt.normalization);

      self.mu = opt.mu;
    end

    w = weight(self, x);
    [a,b] = recurrence(self, n);
  end
  methods(Access=protected)
    p = scale_functions(self, p, n, normalization);
  end

end
