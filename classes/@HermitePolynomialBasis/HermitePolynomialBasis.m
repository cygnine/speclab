classdef HermitePolynomialBasis < OrthogonalPolynomialBasis
  properties(SetAccess=protected)
    mu
  end
  methods
    function self = HermitePolynomialBasis(varargin)
    % self = HermitePolynomialBasis(varargin)
    %
    %     Creates an instance of a Hermite Polynomial spectral basis.

      persistent parser input_parser
      if isempty(parser)
        from labtools import input_parser

        inputs = {'mu', 'domain'};
        defaults = {0, Interval1D([-Inf, Inf])};

        [opt, parser] = input_parser(inputs, defaults, [], varargin{:});
      else
        parser.parse(varargin{:});
        opt = parser.Results;
      end

      opt.standard_domain = Interval1D([-Inf, Inf]);

      self = self@OrthogonalPolynomialBasis(varargin{:});

      self.allowed_function_normalizations{end+1} = ProbabilistFunctionNormalization.instance();
      self.allowed_function_normalizations{end+1} = PhysicistFunctionNormalization.instance();
      self.allowed_function_normalizations{end+1} = ClassicalFunctionNormalization.instance();

      self.mu = opt.mu;
    end

    w = weight(self, x, d);
    [a,b] = recurrence(self, n);
  end
  methods(Access=protected)
    p = scale_functions(self, p, n, normalization);
  end

end
