classdef Chebyshev4PolynomialBasis < JacobiPolynomialBasis
  methods
    function self = Chebyshev4PolynomialBasis(varargin)
    % self = Chebyshev4PolynomialBasis(varargin)
    %
    %     Creates an instance of a Chebyshev4 Polynomial spectral basis.

      persistent all_inputs
      if isempty(all_inputs)
        from labtools import all_inputs
      end
      inputs = {'normalization'};
      defaults = {'normal'};
      opt = all_inputs(inputs, defaults, [], varargin{:});
      opt.alpha = 1/2; 
      opt.beta = -1/2;

      self = self@JacobiPolynomialBasis(opt);
      self.allowed_function_normalizations{end+1} = ClassicalFunctionNormalization.instance();
      self.normalization = self.function_normalization_parser(opt.normalization);
    end
  end

  methods(Access=protected)
    p = scale_functions(self, p, n, normalization);
  end
end
