classdef OrthogonalPolynomialBasis < HilbertBasis
% self = OrthogonalPolynomialBasis({n=0,description='none',
%             fftable=false,domain=Interval1D(),recurrence=@(n) [],
%             normalization='normal', weight_normalization='classical',
%             dim=1, standard_domain=[-1,1]})
%
%     The constructor for an orthogonal polynomial basis set. Is defined
%     almost entirely by the recurrence coefficients. Some custom
%     normalization is possible via optional inputs normalization,
%     weight_normalization, and interval.
%
%     The input "domain" is meant as a user-end input, but the input
%     "standard_domain" is meant to be used in subclassing particular orthogonal
%     polynomial families.
%
%     By default, the input fftable is false, but you may wish to explicitly
%     set it to true; in this case, you are responsible for defining the
%     fft and ifft methods.
  properties
    dim = 1;
    map_to_standard_domain
    map_to_domain
    standard_domain
    domain
  end
  properties(Access=protected)
    default_function_normalization = OrthonormalNormalization.instance();
  end
  properties(SetAccess=protected)
    internal_indexing = ZeroBasedIndexing.instance();
  end
  properties(Access=private)
    recurrence_handle = [];
  end
  methods
    function self = OrthogonalPolynomialBasis(varargin)

      persistent strict_inputs indexing_parser
      if isempty(strict_inputs)
        from labtools import strict_inputs
        from speclab.common import indexing_parser
      end

      %self.default_function_normalization = OrthonormalNormalization.instance();
      %self.default_weight_normalization = ClassicalWeightNormalization.instance();

      inputs = {'recurrence', ...
                'standard_domain', ...
                'indexing', ...
                'normalization', ...
                'weight_normalization', ...
                'dim', ...
                'domain'};
      defaults = {@(n) [], ...
                  Interval1D(), ...
                  ZeroBasedIndexing.instance(), ...
                  'normal', ...
                  'classical', ...
                  1, ...
                  Interval1D()};

      parsed_inputs = strict_inputs(inputs, defaults, [], varargin{:});

      self = self@HilbertBasis(varargin{:});

      % Add allowed normalizations
      self.allowed_function_normalizations{end+1} = MonicNormalization.instance(); 
      self.allowed_function_normalizations{end+1} = OrthonormalNormalization.instance(); 
      self.allowed_weight_normalizations{end+1} = ClassicalWeightNormalization.instance();
      self.allowed_weight_normalizations{end+1} = NaturalWeightNormalization.instance();
      self.allowed_weight_normalizations{end+1} = ProbabilityWeightNormalization.instance();

      self.default_function_normalization = OrthonormalNormalization.instance();
      self.default_weight_normalization = ClassicalWeightNormalization.instance();

      self.recurrence_handle = parsed_inputs.recurrence;

      % Set normalizations
      self.normalization = parsed_inputs.normalization;
      self.weight_normalization = parsed_inputs.weight_normalization;

      % Set up indexing
      self.user_indexing = indexing_parser(parsed_inputs.indexing);
      if isempty(parsed_inputs.internal_indexing)
        self.internal_indexing = indexing_parser(parsed_inputs.indexing);
      end

      % Get domain mapping
      self.standard_domain = parsed_inputs.standard_domain;
      self.domain = Interval1D(parsed_inputs.domain);
    end

    p = evaluate(self,x,n,varargin);
    [x,w] = gauss_quadrature(self,n);
    [x,w] = gauss_radau_quadrature(self,n,varargin);
    [x,w] = gauss_lobatto_quadrature(self,n,varargin);
    [a,b] = recurrence(self, n);
    J = jacobi_matrix(self,N);
    w = weight(self,x);
    k = leading_coefficient(self,n);
    %h = l2_norm(self,n);
    C = monomial_connection(self,N);
    C = inv_monomial_connection(self,N);
    [a,b,c] = mapped_recurrence(self, n);
    lmbda = orthogonal_connection(self,other,N);
    lmbda = self_connection(self,d,N);

    function self = set.domain(self,newdomain)
      self.domain = newdomain;
      self.map_to_standard_domain = self.domain.compute_affine_map(self.standard_domain);
      self.map_to_domain = inv(self.map_to_standard_domain);
    end

    function self = set.user_indexing(self,newindexing)
      if isa(newindexing, 'IndexingRule')
        self.user_indexing = newindexing;
      else
        error('New indexing rule must be object derived from class IndexingRule');
      end
    end

    [n_array, nsize, numeln] = indexing(self,n)
  end

  methods(Access=protected)
    p = scale_functions(self, p, n, normalization)
    w = scale_weight(self,w);
    [x,w] = scale_quadrature(self,x,w);
    p = eval_driver(self, x, a, b, n, d);
  end

end
