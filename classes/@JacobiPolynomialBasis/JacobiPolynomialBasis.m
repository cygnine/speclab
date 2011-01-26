classdef JacobiPolynomialBasis < OrthogonalPolynomialBasis
  properties
    alpha;
    beta;
  end
  methods
    function self = JacobiPolynomialBasis(varargin)
    % self = JacobiPolynomialBasis(varargin)
    %
    %     Creates an instance of a Jacobi Polynomial spectral basis.

      persistent strict_inputs
      if isempty(strict_inputs)
        from labtools import strict_inputs
      end
      inputs = {'alpha', 'beta'};
      defaults = {-0.5, -0.5};
      opt = strict_inputs(inputs, defaults, [], varargin{:});

      self = self@OrthogonalPolynomialBasis(varargin{:});

      [self.alpha, self.beta] = deal(opt.alpha, opt.beta);
      
      if (mod(2*self.alpha,2)==1) && (mod(2*self.beta,2)==1)
        self.fftable = true;
      end
    end

    [a,b] = recurrence(self, n);
    w = weight(self, x);
  end
  methods(Access=protected)
    p = scale_functions(self, p, n, normalization);
  end
end
