classdef HermitePolynomialBasis < OrthogonalPolynomialBasis
  methods
    function self = HermitePolynomialBasis(varargin)
    % self = HermitePolynomialBasis(varargin)
    %
    %     Creates an instance of a Hermite Polynomial spectral basis.

      persistent strict_inputs
      if isempty(strict_inputs)
        from labtools import strict_inputs
      end
      inputs = {'mu'};
      defaults = {0};
      opt = strict_inputs(inputs, defaults, [], varargin{:});

      self = self@OrthogonalPolynomialBasis(varargin{:});

      self.mu = opt.mu;
    end

    w = weight(self, x);
    [a,b] = recurrence(self, n);
  end

end
