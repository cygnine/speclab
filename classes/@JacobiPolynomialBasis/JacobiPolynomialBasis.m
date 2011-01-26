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

      %global packages;
      %inputs = {'alpha', 'beta', 'scale', 'shift', 'dof', 'physical_interval'};
      %defaults = {-0.5, -0.5, 1, 0, 0, false};
      %opt = packages.labtools.input_schema(inputs, defaults, varargin{:});
      %self = basis.WholeBasis(opt);

      %if opt.physical_interval
      %  opt.scale = 1/2*(opt.physical_interval(2) - opt.physical_interval(1));
      %  opt.shift = 1/2*(opt.physical_interval(2) + opt.physical_interval(1));
      %else
      %  opt.physical_interval = [-opt.scale+opt.shift, opt.scale+opt.shift];
      %end
      persistent strict_inputs
      if isempty(strict_inputs)
        from labtools import strict_inputs
      end
      inputs = {'alpha', 'beta'};
      defaults = {-0.5, -0.5};
      opt = strict_inputs(inputs, defaults, [], varargin{:});

      self = self@OrthogonalPolynomialBasis(varargin{:});

      [self.alpha, self.beta] = deal(opt.alpha, opt.beta);
      %self.parameters = params;
      
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
