classdef MultivariateOrthogonalPolynomialBasis < WholeBasis
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
    indexing = @(n) [];
    dimension = 1;
    normalization = [];
    weight_normalization = [];
    map_to_standard_domain
    map_to_domain
    domain
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
    function self = OrthogonalPolynomialBasis(varargin)

      persistent whole_range
      if isempty(whole_range)
        from speclab.common import whole_range
      end

      persistent inparse
      if isempty(inparse)
        inparse = inputParser();
        inparse.KeepUnmatched = true;
        
        inparse.addParamValue('recurrence', @(n) []);
        inparse.addParamValue('standard_domain', Interval1D());
        inparse.addParamValue('normalization', 'normal');
        inparse.addParamValue('weight_recurrence', 'classical');
        inparse.addParamValue('dim', 1); 
        inparse.addParamValue('domain', Interval1D());
      end

      inparse.parse(varargin{:});
      parsed_inputs = inparse.Results;

      %inputs = {'recurrence', 'standard_domain', ...
      %          'normalization', 'weight_normalization', 'dim', 'domain'};
      %defaults = {@(n) [], Interval1D(), 'normal', 'classical', 1, Interval1D()};

      self = self@WholeBasis(varargin{:});

      self.allowed_function_normalizations{end+1} = ClassicalFunctionNormalization.instance(); 
      self.allowed_weight_normalizations{end+1} = ClassicalWeightNormalization.instance();
      self.allowed_weight_normalizations{end+1} = NaturalWeightNormalization.instance();

      self.default_function_normalization = ClassicalFunctionNormalization.instance();
      self.default_weight_normalization = ClassicalWeightNormalization.instance();

      self.normalization = self.function_normalization_parser(parsed_inputs.normalization);
      self.weight_normalization = self.weight_normalization_parser(parsed_inputs.weight_normalization);

      % Get domain mapping
      %self.standard_domain = parsed_inputs.standard_domain;
      %self.domain = Interval1D(parsed_inputs.domain);
      self.indexing = whole_range;
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
