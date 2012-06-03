classdef DiscreteChebyshevPolynomialBasis < OrthogonalPolynomialBasis
  methods
    function self = DiscreteChebyshevPolynomialBasis(varargin)
      persistent parser input_parser struct2cell
      if isempty(parser)
        from labtools import input_parser
        from labtools import struct2fullcell as struct2cell
        inputs = {'N', ...
                  'normalization', ...
                  'standard_domain'};
        defaults = {1, ...
                    'normal', ...
                    []};

        [opt, parser] = input_parser(inputs, defaults, {}, varargin{:});
      else
        parser.parse(varargin{:});
        opt = parser.Results;
      end

      self = self@OrthogonalPolynomialBasis(varargin{:});
      self.standard_domain = Interval1D([0, opt.N-1]);
      self.allowed_function_normalizations{end+1} = ClassicalFunctionNormalization.instance();

      % In case constructor for OrthogonalPolynomialBasis didn't construct it
      self.N = opt.N;
    end

    [a,b] = standard_recurrence(self, n);
    w = weight(self,x);

  end

  methods(Access=protected)
    p = scale_functions(self, p, n, normalization)
  end

end
