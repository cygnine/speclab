classdef JacobiPolynomialBasis < OrthogonalPolynomialBasis
  properties(SetAccess=protected)
    alpha;
    beta;
  end
  methods
    function self = JacobiPolynomialBasis(varargin)
    % self = JacobiPolynomialBasis(varargin)
    %
    %     Creates an instance of a Jacobi Polynomial spectral basis.

      persistent inparse
      if isempty(inparse)
        inparse = inputParser();
        inparse.KeepUnmatched = true;

        inparse.addParamValue('alpha', -0.5);
        inparse.addParamValue('beta', -0.5);
        inparse.addParamValue('normalization', 'normal');
      end

      inparse.parse(varargin{:});
      opt = inparse.Results;

      self = self@OrthogonalPolynomialBasis(varargin{:});
      self.allowed_function_normalizations{end+1} = ClassicalFunctionNormalization.instance();
      self.normalization = opt.normalization;

      [self.alpha, self.beta] = deal(opt.alpha, opt.beta);
      
      if (mod(2*self.alpha,2)==1) && (mod(2*self.beta,2)==1)
        self.fftable = true;
      end

    end

    [a,b] = standard_recurrence(self, n);
    w = weight(self, x);
  end
  methods(Access=protected)
    p = scale_functions(self, p, n, normalization);
  end
end
