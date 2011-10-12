classdef GegenbauerPolynomialBasis < OrthogonalPolynomialBasis
  properties
    lambda
    mu
  end
  properties(Access=private)
    evenbasis
    oddbasis
  end
  methods
    function self = GeneralizedGegenbauerPolynomialBasis(varargin)
    % self = GeneralizedGegenbauerPolynomialBasis(varargin)
    %
    %     Creates an instance of a generalized Gegenbauer Polynomial spectral
    %     basis, which as weight function that is a product of classical weight
    %     functions.

      persistent all_inputs
      if isempty(all_inputs)
        from labtools import all_inputs
      end
      inputs = {'normalization', 'lambda', 'mu'};
      defaults = {'normal', 0, 0};
      opt = all_inputs(inputs, defaults, [], varargin{:});

      opt.alpha = 0-1/2; 
      opt.beta = 0-1/2;

      self = self@OrthogonalPolynomialBasis(opt);
      self.lambda = opt.lambda;
      self.mu = opt.mu;

      self.evenbasis = JacobiPolynomialBasis('alpha', self.lambda-1/2, 'beta', self.mu-1/2, 'normalization', 'classical');
      self.oddbasis = JacobiPolynomialBasis('alpha', self.lambda-1/2, 'beta', self.mu+1/2, 'normalization', 'classical');

      %self.allowed_function_normalizations{end+1} = ClassicalFunctionNormalization.instance();
      self.normalization = self.function_normalization_parser(opt.normalization);
    end
  end

  methods(Access=protected)
    c = constant_cn(self,n,lambda,mu);
    p = scale_functions(self, p, n, normalization);
  end
end
