classdef GegenbauerPolynomialBasis < JacobiPolynomialBasis
  properties
    lambda
  end
  methods
    function self = GegenbauerPolynomialBasis(varargin)
    % self = GegenbauerPolynomialBasis(varargin)
    %
    %     Creates an instance of a Gegenbauer Polynomial spectral basis.

      persistent all_inputs
      if isempty(all_inputs)
        from labtools import all_inputs
      end
      inputs = {'normalization', 'lambda'};
      defaults = {'normal', 0};
      opt = all_inputs(inputs, defaults, [], varargin{:});
      opt.alpha = 0-1/2; 
      opt.beta = 0-1/2;

      self = self@JacobiPolynomialBasis(opt);
      self.lambda = opt.lambda;
      self.allowed_function_normalizations{end+1} = ClassicalFunctionNormalization.instance();
      self.normalization = self.function_normalization_parser(opt.normalization);
    end
  end

  methods(Access=protected)
    p = scale_functions(self, p, n, normalization);
  end
end
