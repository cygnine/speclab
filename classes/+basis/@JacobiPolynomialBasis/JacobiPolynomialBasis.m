global packages;
jac = packages.speclab.OrthgonalPolynomial1D.jacobi;

classdef JacobiPolynomialBasis < basis.WholeBasis
  properties
    parameters;
  end
  methods
    function self = JacobiPolynomialBasis(varargin)
    % self = JacobiPolynomialBasis(varargin)
    %
    %     Creates an instance of a Jacobi Polynomial spectral basis.

      global packages;
      inputs = {'alpha', 'beta', 'scale', 'shift', 'dof', 'physical_interval'};
      defaults = {-0.5, -0.5, 1, 0, 0, false};
      opt = packages.common.input_schema(inputs, defaults, varargin{:});
      self = basis.WholeBasis(opt);

      if opt.physical_interval
        opt.scale = 1/2*(opt.physical_interval(2) - opt.physical_interval(1));
        opt.shift = 1/2*(opt.physical_interval(2) + opt.physical_interval(1));
      else
        opt.physical_interval = [-opt.scale+opt.shift, opt.scale+opt.shift];
      end

      self.interval = opt.physical_interval;
      [params.alpha, params.beta, params.shift, params.scale] = deal(opt.alpha,...
        opt.beta, opt.shift, opt.scale);
      self.parameters = params;
      
      if (mod(2*params.alpha,2)==1) && (mod(2*params.beta,2)==1)
        self.fftable = true;
      end

      self = self.make_vandermonde;
    end

    function output = evaluate(self,x,ns)
      output = self.work_functions.eval.evaluate(x,ns,self.parameters);
    end
    function output = derivative(self,x,ns)
      params = self.parameters;
      params.d = 1;
      output = self.work_functions.eval.evaluate(x,ns,params);
    end

    function value = get.vandermonde(self)
      if isempty(self.vandermonde)
        self = self.make_vandermonde
      end
      value = self.vandermonde;
    end
end
