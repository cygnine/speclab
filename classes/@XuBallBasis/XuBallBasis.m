classdef XuBallBasis < MultivariateOrthogonalPolynomialBasis

  properties
    mu
  end
  properties(Access=private)
    univariate_bases
    recurrence_handle = [];
  end
  methods
    function self = XuBallBasis(varargin)
    % self = XuBallBasis(varargin)
    %
    %     Optional inputs: 'normalization':      normalization for polynomials
    %                      'kappa':              defines weight function
    %                      'mu':                 defines weight function
    %
    %     Creates an instance of a the orthogonal polynomial basis over the unit
    %     ball in 2 dimensions that was defined by Y Xu. The weight function is given by 
    %
    %     w_{mu} = (1 - x^2 - y^2)^{mu-1/2}
    %
    %     for mu > -1/2.
    %
    %     By default this class evaluates the orthonormal polynomials.

      persistent parser input_parser
      if isempty(parser)
        from labtools import input_parser

        inputs = {'mu', ...
                  'indexing', ...
                  'internal_indexing', ...
                  'normalization', ...
                  'weight_normalization'};
        defaults = {1/2, ...
                    [], ...
                    [], ...
                    'normal', ...
                    'classical'};
        [opt, parser] = input_parser(inputs, defaults, [], varargin{:});
      else
        parser.parse(varargin{:});
        opt = parser.Results;
      end
      %inputs = {'normalization', 'kappa', 'mu'};
      %defaults = {'normal', 0, 0};
      %opt = all_inputs(inputs, defaults, [], varargin{:});

      self = self@MultivariateOrthogonalPolynomialBasis(parser.Unmatched);
      self.mu = opt.mu;

      % Right now only generate the basis for y; generate the basis for x as needed
      self.univariate_bases{2} = JacobiPolynomialBasis('alpha', self.mu - 1/2, 'beta', self.mu - 1/2);
      self.univariate_bases{1} = {};

      %self.allowed_function_normalizations{end+1} = ClassicalFunctionNormalization.instance();
      self.allowed_function_normalizations{end+1} = OrthonormalNormalization.instance(); 
      self.allowed_weight_normalizations{end+1} = ProbabilityWeightNormalization.instance();

      self.normalization = self.function_normalization_parser(opt.normalization);
      self.weight_normalization = self.weight_normalization_parser(opt.weight_normalization);

      self.user_indexing = opt.indexing;
      if isempty(opt.indexing)
        user_rules = {};
        user_rules = cell([self.dim 1]);
        user_rules(:) = {ZeroBasedIndexing.instance()};
        self.user_indexing = DirectSumIndexingRule(self.user_indexing, user_rules{:});
      else
        %%%%%%%%%%%%%%%%% For special cases that needs fixing: 
        if isa(self.user_indexing, 'DegreeIndexing');
          self.user_indexing = DegreeIndexing.instance(self.dim);
        end
      end
    end

    p = evaluate(self,x,n);
    h = norm(self,n)
  end

  methods(Access=protected)
    p = scale_functions(self, p, n, normalization);
    w = scale_weight(self, w);
  end

end
