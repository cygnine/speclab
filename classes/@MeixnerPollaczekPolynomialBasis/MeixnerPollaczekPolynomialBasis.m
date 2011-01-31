classdef MeixnerPollaczekPolynomialBasis < OrthogonalPolynomialBasis
  properties
    lambda
    phi
  end
  methods
    function self = MeixnerPollaczekPolynomialBasis(varargin)
    % self = MeixnerPollaczekPolynomialBasis(varargin)
    %
    %     Creates an instance of a Meixner-Pollaczek Polynomial spectral basis.

      persistent all_inputs
      if isempty(all_inputs)
        from labtools import all_inputs
      end
      inputs = {'normalization', 'lambda', 'phi'};
      defaults = {'normal', 1, pi/2};
      opt = all_inputs(inputs, defaults, [], varargin{:});
      opt.alpha = 0-1/2; 
      opt.beta = 0-1/2;

      self = self@OrthogonalPolynomialBasis(opt);
      self.lambda = opt.lambda;
      self.phi = opt.phi;
      self.allowed_function_normalizations{end+1} = ClassicalFunctionNormalization.instance();
      self.normalization = self.function_normalization_parser(opt.normalization);
    end

    [a,b] = recurrence(self, n);
  end

  methods(Access=protected)
    p = scale_functions(self, p, n, normalization);
  end
end
