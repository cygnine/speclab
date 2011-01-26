classdef OrthogonalPolynomialBasis < WholeBasis
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
    indexing = @(n) [];
    dimension = 1;
    normalization = [];
    weight_normalization = [];
    map_to_standard_domain
    map_to_domain
  end
  properties(Access=private)
    recurrence_handle = [];
  end
  methods
    function self = OrthogonalPolynomialBasis(varargin)

      persistent strict_inputs whole_range
      if isempty(strict_inputs)
        from labtools import strict_inputs
        from speclab.common import whole_range
      end

      inputs = {'recurrence', 'standard_domain', ...
                'normalization', 'weight_normalization', 'dim', 'domain'};
      defaults = {@(n) [], Interval1D(), 'normal', 'classical', 1, Interval1D()};

      parsed_inputs = strict_inputs(inputs, defaults, [], varargin{:});

      self = self@WholeBasis(varargin{:});

      self.recurrence_handle = parsed_inputs.recurrence;
      self.normalization = self.function_normalization_parser(parsed_inputs.normalization);
      self.weight_normalization = self.weight_normalization_parser(parsed_inputs.weight_normalization);

      % Get domain mapping
      self.standard_domain = parsed_inputs.standard_domain;
      %self.domain = Orthotope('interval', Interval1D(parsed_inputs.domain));
      self.domain = Interval1D(parsed_inputs.domain);
      self.map_to_standard_domain = self.domain.compute_affine_map(self.standard_domain);
      self.map_to_domain = inv(self.map_to_standard_domain);
      
      self.indexing = whole_range;
    end

    p = evaluate(self,x,n,varargin);
    [x,w] = gauss_quadrature(self,n);
    [x,w] = gauss_radau_quadrature(self,n,varargin);
    [x,w] = gauss_lobatto_quadrature(self,n,varargin);
    [a,b] = recurrence(self, n);
    J = jacobi_matrix(self,N);
    w = weight(self,x);
    k = leading_coefficient(self,n);
    h = l2_norm(self,n);
    C = monomial_connection(self,N);
    C = inv_monomial_connection(self,N);
    [a,b,c] = mapped_recurrence(self, n);
  end

  methods(Access=protected)
    p = scale_functions(self, p, n, normalization)
    w = scale_weight(self,w);
    p = eval_driver(self, x, a, b, n, d);
  end

end