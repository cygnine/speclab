classdef LaguerrePolynomialBasis < OrthogonalPolynomialBasis
  properties
    alpha
  end
  methods
    function self = LaguerrePolynomialBasis(varargin)
    % self = LaguerrePolynomialBasis(varargin)
    %
    %     Creates an instance of a Laguerre Polynomial spectral basis.

      persistent parser input_parser
      if isempty(parser)
        from labtools import input_parser

        inputs = {'alpha', 'domain'};
        defaults = {0, Interval1D([0, Inf])};

        [opt, parser] = input_parser(inputs, defaults, [], varargin{:});
      else
        parser.parse(varargin{:});
        opt = parser.Results;
      end

      opt.standard_domain = Interval1D([0, Inf]);

      self = self@OrthogonalPolynomialBasis(opt);

      self.allowed_function_normalizations{end+1} = ClassicalFunctionNormalization.instance();

      self.alpha = opt.alpha;
    end

    w = weight(self, x);
    [a,b] = recurrence(self, n);
  end
  methods(Access=protected)
    p = scale_functions(self, p, n, normalization);
  end

end
