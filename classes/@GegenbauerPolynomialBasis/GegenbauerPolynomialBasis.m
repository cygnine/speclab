classdef GegenbauerPolynomialBasis < JacobiPolynomialBasis
  properties
    lambda
  end
  methods
    function self = GegenbauerPolynomialBasis(varargin)
    % self = GegenbauerPolynomialBasis(varargin)
    %
    %     Creates an instance of a Gegenbauer Polynomial spectral basis.

      persistent parser input_parser
      if isempty(parser)
        from labtools import input_parser

        inputs = {'normalization', 'lambda'};
        defaults = {'classical', 0};
        [results, parser] = input_parser(inputs, defaults, [], varargin{:});
      else
        parser.parse(varargin{:});
        opt = parser.Results;
      end

      self = self@JacobiPolynomialBasis('alpha', opt.lambda - 1/2, 'beta', opt.lambda - 1/2);
      self.lambda = opt.lambda;
      self.allowed_function_normalizations{end+1} = ClassicalFunctionNormalization.instance();
      self.normalization = self.function_normalization_parser(opt.normalization);
    end
  end

  methods(Access=protected)
    p = scale_functions(self, p, n, normalization);
  end
end
