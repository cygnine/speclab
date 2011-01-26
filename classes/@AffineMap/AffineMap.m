classdef AffineMap
  properties 
    domain = [];
    range = [];
    A = [];
    b = [];
  end
  properties(SetAccess=private)
    m = 0;
    n = 0;  
  end
  properties(Access=private)
    Aaug = 0;
  end
  methods 
    function self = AffineMap(A, b, varargin)
    % AffineMap -- A Euclidean affine map object
    %
    % self = AffineMap(A, b, {domain=[], range=[]})
    %
    %     Returns an instance of an affine map, that is a map f: x ---> A*x + b,
    %     where x \in R^m, and f(x) \in R^n. A is an (n x m) matrix, and b is an
    %     (n x 1) vector. Two methods of specifying the map are supported:
    %
    %     1.) map = AffineMap(A, b)
    %         Explicitly defines the map coefficients. Note that if A is not
    %         invertible then the inverse map cannot be computed/evaluated.
    %     2.) map = AffineMap([], [], 'domain', domain, 'range', range)
    %         The input variables domain and range are Orthotope object-type
    %         specifications. In this case, the domain and range must have the
    %         same dimension. Note that if either the domain or range contain
    %         (semi-)infinite slices, then there is not a unique affine map
    %         connecting the domain and the range -- in this case a scaling of 1
    %         is chosen, and an appropriate shift is computed. To override this
    %         behavior, use calling syntax 1.).
      
      persistent strict_inputs 
      if isempty(strict_inputs) 
        from labtools import strict_inputs 
      end

      opt = strict_inputs({'domain', 'range'}, {[], []}, [], varargin{:});

      if isempty(A)
        % Then construct affine map using given orthotopes.
        if not(isa(opt.domain, 'Orthotope'))
          if isempty(opt.domain)
            opt.domain = Orthotope();
          else
            opt.domain = Orthotope(opt.domain);
          end
        end
        if not(isa(opt.range, 'Orthotope'))
          if isempty(opt.range)
            opt.range = Orthotope();
          else
            opt.range = Orthotope(opt.range);
          end
        end
        if opt.domain.dimension ~= opt.range.dimension
          error('Affine map cannot be specified with domain and range of different dimensions');
        end
        self.m = self.domain.dimension;
        self.n = self.range.dimension;
        self.A = eye(self.m);
        self.b = zeros([self.n 1]);
        for q = 1:self.n
          self.b(q) = opt.range.slices{q}.centroid - opt.domain.slices{q}.centroid;
        end
      else
        % Pretty straightforward
        [self.n, self.m] = size(A);
        if length(b) ~= self.n
          error('Specified coefficients don''t have conforming size');
        end
        self.A = A;
        self.b = b;
      end

      %self.domain = Orthotope('dimension', self.m, 'interval', [-Inf, Inf]);
      %self.range = Orthotope('dimension', self.n, 'interval', [-Inf, Inf]);
      self.Aaug = [self.A self.b; zeros([1 self.m]) 1];
    end

    c = evaluate(self, x);
    b = compose(self, a);
    m = inv(self);
  end
end
