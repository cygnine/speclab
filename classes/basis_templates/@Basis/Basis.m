classdef Basis
% self = basis({dof=0,description='',fftable=false, interval=[0,0]})
%
%     The generic basis class for all spectral expansions. This is a parent
%     class used by other classes. It is not meant to be useful as a
%     stand-alone class.
  properties
    description = [];
    fftable = false;
    normalization = [];
  end
  properties(Abstract=true)
    domain
    standard_domain
  %end
  %properties(Abstract=true,SetAccess=protected)
    user_indexing
  end
  properties(SetAccess=protected)
    internal_indexing 
    internal_indexset
  end
  properties(Access=protected)
    allowed_function_normalizations = {};
    default_function_normalization
  end
  methods(Abstract=true)
    indexing
  end
  methods
    function self = Basis(varargin)
      persistent strict_inputs
      if isempty(strict_inputs)
        from labtools import strict_inputs
      end

      inputs = {'description','fftable','domain','standard_domain'};
      defaults = {[], false, Interval1D(),Interval1D()};
      parsed_inputs = strict_inputs(inputs, defaults, [], varargin{:});

      self.description = parsed_inputs.description;
      self.fftable = parsed_inputs.fftable;
    end
    function[self] = set.internal_indexing(self,inp)
      self.internal_indexset = inp.image;
      self.internal_indexing = inp;
    end
    function[self] = set.normalization(self, inp)
      self.normalization = self.function_normalization_parser(inp);
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

    %inds = natural_indexing(self, subs)

  end
  methods(Access=protected)
    obj = function_normalization_parser(self, inp);
    obj = indexing_parser(self, inp);
  end
end
