classdef LegendrePolynomialBasis < JacobiPolynomialBasis
% self = LegendrePolynomialBasis()
%
%     Constructor for the LegendrePolynomialBasis class. The default is an
%     orthogonal polynomial family over [-1,1] with weight function w(x) = 1,
%     orthonormal normalization. All customizations given by
%     OrthogonalPolynomialBasis are supported. New normalization types supported
%     by this subclass are:
%
%        'classical' -- p_n(+1) = 1, on [-1,1]. For domain [a,b], p_n(b) = 1.
  methods
    function self = LegendrePolynomialBasis(varargin)
    % self = LegendrePolynomialBasis(varargin)
    %
    %     Creates an instance of a Legendre Polynomial spectral basis.

      persistent all_inputs
      if isempty(all_inputs)
        from labtools import all_inputs
      end
      inputs = {'normalization'};
      defaults = {'normal'};
      opt = all_inputs(inputs, defaults, [], varargin{:});

      opt.alpha = 0; 
      opt.beta = 0;

      self = self@JacobiPolynomialBasis(opt);
      self.normalization = opt.normalization;
    end
  end

  methods(Access=protected)
    p = scale_functions(self, p, n, normalization);
  end
end
