classdef LegendrePolynomialBasis < JacobiPolynomialBasis
  methods
    function self = LegendrePolynomialBasis(varargin)
    % self = LegendrePolynomialBasis(varargin)
    %
    %     Creates an instance of a Legendre Polynomial spectral basis.

      persistent all_inputs
      if isempty(all_inputs)
        from labtools import all_inputs
      end
      inputs = {};
      defaults = {};
      opt = all_inputs(inputs, defaults, [], varargin{:});
      opt.alpha = 0; 
      opt.beta = 0;

      self = self@JacobiPolynomialBasis(opt);
    end
  end

  methods(Access=protected)
    p = scale_functions(self, p, n, normalization);
  end
end
