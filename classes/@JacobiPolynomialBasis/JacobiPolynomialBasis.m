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

      persistent input_parser parser
      if isempty(parser)
        from labtools import input_parser

        inputs = {'alpha', ...
                  'beta', ...
                  'normalization'};
        defaults = {-0.5, ...
                    -0.5, ...
                    'normal'};

        [opt, parser] = input_parser(inputs, defaults, {}, varargin{:});

        %inparse = inputParser();
        %inparse.KeepUnmatched = true;

        %inparse.addParamValue('alpha', -0.5);
        %inparse.addParamValue('beta', -0.5);
        %inparse.addParamValue('normalization', 'normal');
      else
        parser.parse(varargin{:});
        opt = parser.Results;
      end

      self = self@OrthogonalPolynomialBasis(varargin{:});
      self.allowed_function_normalizations{end+1} = ClassicalFunctionNormalization.instance();

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
