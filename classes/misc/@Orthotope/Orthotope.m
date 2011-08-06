classdef Orthotope
% Orthotope({boundaries={}, dimension=1, interval=Interval1D('interval', [-1,1])})
%
%     Creates a tensor-product hypercube from one-dimensional intervals. If
%     boundaries is non-empty: boundaries should be a cell array, where each
%     entry in the cell is an Interal1D-like object. (I.e. either an Interval1D
%     object, or a two-element vector.) The dimension of the Orthotope is then
%     set to be length(boundaries), and each dimension uses the interval
%     specified in the cell array boundaries.
%
%     If boundaries is empty: creates a hypercube in 'dimension'-dimensional
%     space using the input 'interval' as the template.
%
% Orthotope Properties:
%   dimension - (defining) The dimension of the orthotope
%   slices - (defining) Cell array containing Interval1D objects
% 
% Orthotope Methods:
%   compute_affine_map - Computes affine map between this orthoope and another

  properties
    dimension = 1;
    slices = [];
  end
  properties(Access=private)
    default_interval = Interval1D([-1,1]);
  end
  methods
    function self = Orthotope(varargin)
      persistent inparse
      if isempty(inparse)
        inparse = inputParser();
        inparse.KeepUnmatched = false;

        inparse.addParamValue('boundaries', {});
        inparse.addParamValue('dimension',1);
        inparse.addParamValue('interval', self.default_interval);
      end

      %inputs = {'boundaries', 'dimension', 'interval'};
      %defaults = {{}, 1, self.default_interval};;

      inparse.parse(varargin{:});
      temp = inparse.Results;

      if isempty(temp.boundaries)
        % This just tensor-product one-dimensional interval
        if not(isa(temp.interval, 'Interval1D'))
          temp.interval = Interval1D(temp.interval);
        end
        for q = 1:temp.dimension
          self.slices{q} = temp.interval;
        end
      else
        % boundaries has something in it...time to parse
        for q = 1:length(temp.boundaries)
          if isa(temp.boundaries{q}, 'Interval1D')
            self.slices{q} = temp.boundaries{q};
          else
            self.slices{q} = Interval1D(temp.boundaries{q});
          end
        end
      end
      self.dimension = length(self.slices);
    end
    
    m = compute_affine_map(self, other);
  end
end
