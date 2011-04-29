classdef Basis
% self = basis({dof=0,description='',fftable=false, interval=[0,0]})
%
%     The generic basis class for all spectral expansions. This is a parent
%     class used by other classes. It is not meant to be useful as a
%     stand-alone class.
  properties
    description = [];
    fftable = false;
  end
  properties(Abstract=true)
    domain
    standard_domain
    indexing 
  end
  properties(Abstract=true,Access=protected)
    default_indexing
  end
  properties(Access=protected)
    allowed_function_normalizations = {};
    allowed_weight_normalizations = {};
  %end
  %properties(Access=protected)
    default_function_normalization
    default_weight_normalization
  end
  methods
    function self = Basis(varargin)
      persistent strict_inputs
      if isempty(strict_inputs)
        from labtools import strict_inputs
      end

      inputs = {'description','fftable','domain','standard_domain'};
      defaults = {[], false, Interval1D(),Interval1D()};
      opt = strict_inputs(inputs, defaults, [], varargin{:});

      self.description = opt.description;
      self.fftable = opt.fftable;
    end

    function output = range(self,N)
    % range -- Returns indices for the first N basis elements
    %
    % output = range(self,N)
    %
    %     Given an integer N, this function returns indices for the first N
    %     basis elements.
      output = self.indexing.from_naturals(1:N);
    end

    obj = function_normalization_parser(self, inp);
    obj = weight_normalization_parser(self, inp);
    %inds = natural_indexing(self, subs)

  end
end
