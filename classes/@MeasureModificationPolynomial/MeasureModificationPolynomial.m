classdef MeasureModificationPolynomial < FactoredPolynomial
% MeasureModificationPolynomial -- A polynomial for use in (positive) 
%                                  measure modifications
%
% self = MeasureModificationPolynomial(interval, roots, 
%                                      [[ leading_coefficient = 1 ]] )
%
%     Creates a FactoredPolynomial-derived instance that is designed to be
%     semipositive on 'interval'. 'interval' is either (a) a length-2
%     numeric array, or (b) an Interval1D instance that defines a univariate
%     interval. Upon verification of the interval, the constructor ensures
%     that the input roots define a polynomial that is semipositive on the
%     interval.  A tolerance is defined (1e-10) to make this happen (e.g.
%     for quadratic roots, or complex-conjugate roots).
%
% self = MeasureModificationPolynomial(interval, exterior_roots,
%                                      conjugate_roots, quadratic_roots, 
%                                      [[ leading_coefficient = 1 ]])
%
%     Explicitly delineates the types of roots: 
%        exterior_roots -- real-valued roots that lie outside interval
%        conjugate_roots -- complex-valued roots whose c.c. root is added
%        quadratic_roots -- real-valued roots inside the interval that are 
%                           quadratic roots

% Matlab has perhaps the stupidest rule possible about excluding conditional
% calls to superclass constructors. Hence the following ridiculous
% convolution of logical statements.
  properties
    exterior_roots
    conjugate_roots
    quadratic_roots 
    interval
    tol = 1e-10;
  end
  properties(Access=protected)
    num_exterior = 0;
    num_conjugate = 0;
    num_quadratic = 0;
  end
  methods
    function[self] = MeasureModificationPolynomial(varargin)
      if nargin < 4
        myargs = varargin(2:end);
      else
        myargs = {[]};
      end

      self = self@FactoredPolynomial(myargs{:});

      if nargin < 4  % Calling syntax is FactoredPolynomial type
        self.interval = Interval1D(varargin{1});
      else
        self.interval = Interval1D(varargin{1});
        if nargin < 5
          varargin{5} = 1;
        end
        % Given factorization explicitly
        self.exterior_roots = varargin{2};
        self.conjugate_roots = varargin{3};
        self.quadratic_roots = varargin{4};
        self.leading_coefficient = varargin{5};
        self = self.fullroots_from_modified_roots();
      end

      self = self.verify_positivity();
    end

    function[self] = set.exterior_roots(self, newroots)
      % Update num_exterior
      self.exterior_roots = newroots(:);
      self.num_exterior = numel(self.exterior_roots);
    end

    function[self] = set.conjugate_roots(self, newroots)
      % Update num_conjugate
      self.conjugate_roots = newroots(:);
      self.num_conjugate = numel(self.conjugate_roots);
    end

    function[self] = set.quadratic_roots(self, newroots)
      % Update num_quadratic
      self.quadratic_roots = newroots(:);
      self.num_quadratic = numel(self.quadratic_roots);
    end

    varargout = plot(self, varargin)
  end
  methods(Access=protected)
    self = fullroots_from_modified_roots(self);
    self = verify_positivity(self,newroots);
  end
end
