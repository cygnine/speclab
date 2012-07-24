classdef GeneralizedGegenbauerPolynomialBasis < OrthogonalPolynomialBasis
  properties
    lambda
    mu
  end
  properties(Access=private)
    evenbasis
    oddbasis
  end
  methods
    function self = GeneralizedGegenbauerPolynomialBasis(varargin)
    % self = GeneralizedGegenbauerPolynomialBasis(varargin)
    %
    %     Creates an instance of a generalized Gegenbauer Polynomial spectral
    %     basis, which as weight function that is a product of classical weight
    %     functions.

      persistent input_parser parser
      if isempty(parser)
        from labtools import input_parser

        inputs = {'lambda', ...
                  'mu', ...
                  'normalization'};
        defaults = {0, ...
                    0, ...
                    'normal'};

        [opt, parser] = input_parser(inputs, defaults, {}, varargin{:});

      else
        parser.parse(varargin{:});
        opt = parser.Results;
      end

      self = self@OrthogonalPolynomialBasis(opt);
      self.lambda = opt.lambda;
      self.mu = opt.mu;
    end

    [a,b] = standard_recurrence(self, n);
    w = weight(self, x);
  end

  methods(Access=protected)
    p = scale_functions(self, p, n, normalization);
  end
end
