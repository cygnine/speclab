classdef Interval1D
% Interval1D(interval {property-key-value-pairs})
% 
%     A one-dimensional interval object. Mainly used for computing affine maps.
%     Is not classified as either open, closed, or half-open. 
%
%     This class is implemented to always have an 'anchor' affine map connecting
%     [-1,1] to the specified interval. If the interval has finite-length, this
%     is well-defined. For infinite-length intervals, additional information
%     must be supplied.
%
%     The optional inputs 'centroid' and 'scale' are ignored if the interval has
%     finite length. If the interval has infinite length, then 'centroid'
%     corresponds to the 'origin' location of the interval, and 'scale' is the
%     distance between the origin and the image of 1. These quantities are only
%     necessary if one wishes to form an affine map between the interval and
%     some other interval.
%
%     E.g. : Interval1D([-Inf, Inf], 'centroid', 3, 'scale', 0.5)
%              represents the whole real line, where x=3 is considered the
%              origin, and [2.5, 3.5] is considered the image of [-1,1] under an
%              affine map.
%
% Interval1D Properties:
%   interval - (assignable) 
%   centroid - (assignable) The interval 'center' (default: 0)
%   scale - (assignable) The distance from the image of 0 to the image of 1 under an affine map that defines the interval
%   length - The length of the interval
%   map_to_standard_interval - An affine map between the standard interval [-1,1] and this interval
%
% Interval1D Methods:
%   linspace - equispaced points on the interval
%   compute_affine_mape - Computes affine map between this interval and another

  properties(SetAccess=private)
    centroid = 0; 
    scale = 1;
    length = 0;
    interval = [0, 0];
  end
  properties(SetAccess=private)
    map_to_standard_interval = AffineMap(1, 0);
  end
  methods
    function self = Interval1D(interval, varargin)
      persistent inparse 
      if isempty(inparse)
        inparse = inputParser();
        inparse.KeepUnmatched = false;

        inparse.addParamValue('centroid', []);
        inparse.addParamValue('scale',1);
      end

      % Gives default value of 'interval' in case constructor is called with no
      % inputs
      if not(exist('interval'))
        interval = [-1, 1];
      end

      % If an Interval1D instance is given as input, just throw it back.
      if isa(interval, 'Interval1D');
        self = interval;
        return
      end

      inparse.parse(varargin{:});
      temp = inparse.Results;

      self.interval = interval(:).';
      self.centroid = temp.centroid;

      if length(self.interval) ~= 2
        error('Input must be a two-element vector.');
      end
      % If user tries to create e.g. interval [+1, -1], reverse order
      if diff(self.interval)<0
        fprintf(['Warning: specified interval has negative length -- I am \\' ...
                 'reversing the order of the defined interval.\n']);
        temp = self.interval(1);
        self.interval(1) = self.interval(2); 
        self.interval(2) = temp;
      end

      self.length = diff(self.interval);
      self.centroid = mean(self.interval);
      self.scale = self.length/2;

      % If the interval has infinite length
      if any(isinf(self.interval))
        % Use the given centroid value, if any
        self.scale = temp.scale;
        self.centroid = temp.centroid;
        if isinf(self.interval(1)) & isinf(self.interval(2))
          % The interval is [-Inf, Inf]
          if isempty(temp.centroid)
            self.centroid = 0;
          end
        elseif isinf(self.interval(1))
          % The interval is [-Inf, something]
          if isempty(temp.centroid)
            self.centroid = self.interval(2) - self.scale;
          end
        elseif isinf(self.interval(2))
          % The interval is [something, Inf]
          if isempty(temp.centroid)
            self.centroid = self.interval(1) + self.scale;
          end
        end
      end

      self.map_to_standard_interval = AffineMap(1/self.scale, -1/self.scale*self.centroid);
    end

    x = linspace(self, N);
    m = compute_affine_map(self, varargin)
  end
end
