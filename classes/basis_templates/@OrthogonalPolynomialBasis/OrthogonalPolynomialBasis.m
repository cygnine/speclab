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
%
% OrthogonalPolynomialBasis Properties:
%   indexing - IndexingRule-derived object specifying a map from user-end indexing to the natural numbers 1, 2, ...
%   normalization - Normalization-derived object specifying the function normalizations
%   domain - User-end domain after affine map from standard_domain
%   standard_domain - Standard (internal use) univariate domain
%   map_to_domain - The affine map from standard_domain to domain
%   map_to_standard_domain - The affine map from domain to standard_domain
%
% OrthogonalPolynomialBasis Methods:
%   evaluate - Evaluates polynomials
%   gauss_quadrature - Gauss quadrature rule
%   gauss_radau_quadrature - Gauss-Radau quadrature rule
%   gauss_lobatto_quadrature - Gauss-Lobatto quadrature rule
%   recurrence - Recurrence coefficients
%   jacobi_matrix - Jacobi matrix
%   weight - Evaluates the weight function 
%   norm - Evaluates the l^2 norm of the polynomials
%   leading_coefficient - Returns leading coefficients for polynomials
%   monomial_connection - Returns the connection matrix between this family and the monomials
%   inv_monomial_connection - Returns the connection matrix between the monomials and this family
%   mapped_recurrence - Recurrence coefficients for the polynomials after mapping to self.domain
%   orthogonal_connection - Returns the connection matrix between this family and another OrthogonalPolynomialBasis family
%   derivative_expansion - expansion coefficients of polynomial derivatives
  properties
    dim = 1; % Make this hidden
    map_to_standard_domain % A map from domain to standard_domain derived from domain and standard_domain
    map_to_domain % A map from standard_domain to domain derived from domain and standard_domain
    standard_domain % { Interval1D([-1, 1]) } | Other_Interval1D_object - Standard univariate interval for basis
  end
  properties(Access=private)
    recurrence_handle = [];
  end
  methods
    function self = OrthogonalPolynomialBasis(varargin)

      persistent inparse
      if isempty(inparse)
        inparse = inputParser();
        inparse.KeepUnmatched = true;

        inparse.addParamValue('recurrence', @(n) []);
        inparse.addParamValue('standard_domain', Interval1D());
        inparse.addParamValue('indexing', ZeroBasedIndexing.instance());
        inparse.addParamValue('internal_indexing', []);
        inparse.addParamValue('normalization', 'normal');
        inparse.addParamValue('weight_normalization', 'classical');
        inparse.addParamValue('dim', 1);
        inparse.addParamValue('domain', Interval1D());
      end

      self = self@HilbertBasis(varargin{:});

      % Defaults for OrthogonalPolynomialBasis:
      self.internal_indexing = ZeroBasedIndexing.instance();
      self.default_indexing_rule = ZeroBasedIndexing.instance();
      self.default_weight_normalization = ClassicalWeightNormalization.instance();

      % Prescribe allowed normalizations
      self.allowed_function_normalizations{end+1} = MonicNormalization.instance(); 
      self.allowed_function_normalizations{end+1} = OrthonormalNormalization.instance(); 
      self.allowed_weight_normalizations{end+1} = ClassicalWeightNormalization.instance();
      self.allowed_weight_normalizations{end+1} = NaturalWeightNormalization.instance();
      self.allowed_weight_normalizations{end+1} = ProbabilityWeightNormalization.instance();

      % Parse inputs
      inparse.parse(varargin{:});
      parsed_inputs = inparse.Results;

      % Deal with given inputs
      self.recurrence_handle = parsed_inputs.recurrence;

      % Set normalizations
      self.normalization = parsed_inputs.normalization;
      self.weight_normalization = parsed_inputs.weight_normalization;

      % Set up indexing
      self.user_indexing = parsed_inputs.indexing;
      self.internal_indexing = parsed_inputs.internal_indexing;

      % Get domain mapping
      self.standard_domain = Interval1D(parsed_inputs.standard_domain);
      self.domain = Interval1D(parsed_inputs.domain);
    end

    p = evaluate(self,x,n,varargin);
    [x,w] = gauss_quadrature(self,n);
    [x,w] = gauss_radau_quadrature(self,n,varargin);
    [x,w] = gauss_lobatto_quadrature(self,n,varargin);
    [a,b] = recurrence(self, n);
    J = jacobi_matrix(self,N);
    w = weight(self,x);
    h = norm(self,n)
    k = leading_coefficient(self,n);
    C = monomial_connection(self,N);
    C = inv_monomial_connection(self,N);
    [a,b,c] = mapped_recurrence(self, n);
    lmbda = orthogonal_connection(self,other,N);
    lmbda = self_connection(self,d,N);
    lmbda = derivative_expansion(self, N, d);
  end

  methods(Access=protected)
    p = scale_functions(self, p, n, normalization)
    w = scale_weight(self,w);
    [x,w] = scale_quadrature(self,x,w);
    %p = eval_driver(self, x, a, b, n, d);
  end

  methods(Static)
    [x,w] = gauss_quadrature_driver(a,b)
    p = evaluate_driver(x, alpha, beta, n, d)
  end

end
