classdef Interval1D
% Interval1D({interval=[], centroid=[]})
% 
%     A one-dimensional interval object. Mainly used for computing affine maps.
%     Is not classified as either open, closed, or half-open. The `centroid'
%     input is ignored unless the given interval has length infinity. Note that
%     in this case, the `centroid' is meaningless -- it is only used to
%     construct affine maps.
  properties
    length = 0;
    interval = [0, 0];
    centroid = 0;
  end
  methods
    function self = Interval1D(varargin)
      persistent strict_inputs
      if isempty(strict_inputs)
        from labtools import strict_inputs
      end

      inputs = {'interval', 'centroid'};
      defaults = {[0,0], []};
      % All I want to do is:
      % self = strict_inputs(inputs, defaults, [], varargin{:});
      % Are you kidding me, Matlab?
      temp = strict_inputs(inputs, defaults, [], varargin{:});
      tempnames = fieldnames(temp);
      for q = 1:length(tempnames)
        self = setfield(self, tempnames{q}, getfield(temp, tempnames{q}));
      end
      self.interval = self.interval(:).';  % Force row vector

      if length(self.interval) ~= 2
        error('Input must be a two-element vector.');
      end
      if diff(self.interval)<0
        fprintf(['Warning: specified interval has negative length -- I am \\' ...
                 'reversing the order of the defined interval.\n']);
        temp = self.interval(1);
        self.interval(1) = self.interval(2); 
        self.interval(2) = temp;
      end
      self.length = diff(self.interval);
      self.centroid = mean(self.interval);
      if any(isinf(self.interval))
        self.centroid = temp.centroid;
        if isinf(self.interval(1)) & isinf(self.interval(2))
          % The interval is [-Inf, Inf]
          if isempty(temp.centroid)
            self.centroid = 0;
          end
        elseif isinf(self.interval(1))
          % The interval is [-Inf, something]
          if isempty(temp.centroid)
            self.centroid = self.interval(2) - 1;
          end
        elseif isinf(self.interval(2))
          % The interval is [something, Inf]
          if isempty(temp.centroid)
            self.centroid = self.interval(1) + 1;
          end
        end
      end
    end

  end
end
