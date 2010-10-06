classdef Hypercube
% Hypercube({boundaries={}, dim=1, interval=Interval1D('interval', [-1,1])})
%
%     Creates a tensor-product hypercube from one-dimensional intervals. If
%     boundaries is non-empty: boundaries should be a cell array, where each
%     entry in the cell is an Interal1D-like object. (I.e. either an Interval1D
%     object, or a two-element vector.) The dimension of the Hypercube is then
%     set to be length(boundaries), and each dimension uses the interval
%     specified in the cell array boundaries.
%
%     If boundaries is empty: creates a hypercube in dim-dimensional space using
%     the third input interval as the template.

  properties
    dim = 1;
    boundaries = [];
  end
  methods
    function self = Hypercube(varargin)
      persistent strict_inputs default_interval
      if isempty(strict_inputs)
        from labtools import strict_inputs
        default_interval = Interval1D('interval', [-1, 1]);
      end

      inputs = {'boundaries', 'dim', 'interval'};
      defaults = {{}, 1, default_interval};
      temp = strict_inputs(inputs, defaults, [], varargin{:});

      if isempty(temp.boundaries)
        % This just tensor-product one-dimensional interval
        if not(isa(temp.interval, 'Interval1D'))
          temp.interval = Interval1D('interval', temp.interval);
        end
        for q = 1:temp.dim
          self.boundaries{q} = temp.interval;
        end
      else
        % boundaries has something in it...time to parse
        for q = 1:length(temp.boundaries)
          if isa(temp.boundaries{q}, 'Interval1D')
            self.boundaries{q} = temp.boundaries{q};
          else
            self.boundaries{q} = Interval1D('interval', temp.boundaries{q});
          end
        end
      end
      dim = length(self.boundaries);
    end
  end
end
