classdef MultivariateOrthogonalPolynomialBasis < HilbertBasis
% self = MultivariateOrthogonalPolynomialBasis({n=0,description='none',
%             domain=Interval1D(),
%             normalization='normal', weight_normalization='classical',
%             dim=1, affine_map})
%
%     The constructor for a multivariate orthogonal polynomial basis set. This
%     is mainly a base class from which derivations specify the precise nature
%     of the basis. 
%
%     Usually the standard domain is defined by the derived classes, but the
%     input 'affine_map' is a user-end method to map the standard domains to
%     other regions.

  properties
    dim = 2;
    map_to_standard_domain
    map_to_domain
    standard_domain
  end
  properties(Access=private)
    recurrence_handle = [];
  end
  properties(Access=protected)
    %allowed_function_normalizations;
    %allowed_weight_normalizations;
    %default_function_normalization = OrthonormalNormalization.instance();
    %default_weight_normalization = ClassicalWeightNormalization.instance();
  end
  methods
    function self = MultivariateOrthogonalPolynomialBasis(varargin)

      persistent parser input_parser
      if isempty(parser)
        from labtools import input_parser
        inputs = {'recurrence', ...
                  'standard_domain', ...
                  'normalization', ...
                  'weight_normalization', ...
                  'domain', ...
                  'indexing', ...
                  'internal_indexing'};
        defaults = {@(n) [], ...
                    Interval1D(), ...
                    'classical', ...
                    'classical', ...
                    Interval1D(),...
                    [], ...
                    []};
        [parsed_inputs, parser] = input_parser(inputs, defaults, [], varargin{:});
      else
        parser.parse(varargin{:});
        parsed_inputs = parser.Results;
      end

      %self = self@Basis(varargin{:});
      self = self@HilbertBasis(varargin{:});

      % Get domain mapping
      self.standard_domain = Interval1D(parsed_inputs.standard_domain);
      self.domain = Interval1D(parsed_inputs.domain);

      self.allowed_function_normalizations{end+1} = ClassicalFunctionNormalization.instance(); 
      self.allowed_weight_normalizations{end+1} = ClassicalWeightNormalization.instance();
      self.allowed_weight_normalizations{end+1} = NaturalWeightNormalization.instance();

      self.default_function_normalization = ClassicalFunctionNormalization.instance();
      self.default_weight_normalization = ClassicalWeightNormalization.instance();

      self.normalization = self.function_normalization_parser(parsed_inputs.normalization);
      self.weight_normalization = self.weight_normalization_parser(parsed_inputs.weight_normalization);

      % Set up indexing
      self.default_indexing_rule = TotalGRevLexIndexing.instance(2);

      self.user_indexing = parsed_inputs.indexing;
      self.internal_indexing = parsed_inputs.internal_indexing;

    end

    p = evaluate(self,x,n,varargin);
    w = weight(self,x);
    h = l2_norm(self,n);
  end

  methods(Access=protected)
    p = scale_functions(self, p, n, normalization)
    w = scale_weight(self,w);
  end

end
