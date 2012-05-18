classdef EuclideanBall
% EuclideanBall([[dim=1, center=0, radius=1]])
% 
%     A dim-dimensional Euclidean ball. Mainly used for computing affine maps.
%     Is not classified as either open or closed.
%
%     This class is implemented to always have an 'anchor' affine map connecting
%     the unit ball to the specified ball. If the ball has finite radius, this
%     is well-defined. For infinite-sized balls, additional information
%     must be supplied. Examples:
%
%         EuclideanBall('dim', 3)   % 3-dimensional origin-centered unit ball
%         EuclideanBall('dim', 4, 'center', [1 2 3 1]) % 4-dimensional unit-radius ball centered at (1,2,3,1)
%         EuclideanBall('dim', 2, 'radius', Inf) % All of R^2, the affine map is just the identity
%         EuclideanBall('dim', 2, 'radius', Inf, 'center', [-2,3], 'scale', 2) % All of R^2, the affine map from the unit ball maps to the unit ball centered at (-2,3) with radius 2.
%         EuclideanBall('dim', 2, 'scale', 1) % The implied radius is 1, which overrides the scale input
%
% EuclideanBall Properties:
%   dimension - dimension of Euclidean space
%   center - The ball center (default: origin)
%   radius - The ball radius
%   scale - The distance from the image of 0 to the image of 1 under an affine map that defines the unit ball. Same as 'radius' for finite-radius balls.
%   map_to_standard_domain - An affine map between the standard domain (the unit ball) and this domain
%
% Interval1D Methods:
%   compute_affine_map - Computes affine map between this ball and another

  properties(SetAccess=private)
    dimension = 1; 
    center = 0;
    radius = 1;
    scale = 1;
    map_to_standard_domain = AffineMap(1, 0);
  end
  methods
    function self = EuclideanBall(varargin)
      persistent parser input_parser 
      if isempty(parser)
        from labtools import input_parser
        inputs = {'dim', ...
                  'center', ...
                  'radius', ...
                  'scale'};
        defaults = {1, ...
                    0, ...
                    1, ...
                    1};
        validators = {@(dim) assert(length(dim(:))==1, 'Dimension must be a scalar'), ...
                      [], ...
                      @(radius) assert((length(radius(:))==1) && (radius > 0), 'Radius must be a positive scalar'), ...
                      @(scale) assert((length(scale(:))==1) && (scale > 0), 'Scale parameter must be a positive scalar')};
        [opt, parser] = input_parser(inputs, defaults, validators, varargin{:});
      else
        parser.parse(varargin{:});
        opt = parser.Results;
      end

      % Deal user inputs into instance properties
      self.dimension = opt.dim;
      self.radius = opt.radius;
      if opt.center==0
        opt.center = zeros([1 self.dimension]);
      end
      assert(length(opt.center(:))==self.dimension, ['Center of ball must be a ' num2str(self.dimension) '-sized vector']);
      self.center = opt.center;
      if isinf(self.radius)
        self.scale = opt.scale;
      else
        self.scale = self.radius;
      end

      self.map_to_standard_domain = AffineMap(1/self.scale*eye(self.dimension), -1/self.scale*self.center);
    end

    m = compute_affine_map(self, varargin)
  end
end
