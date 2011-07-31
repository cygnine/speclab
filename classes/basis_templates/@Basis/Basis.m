classdef Basis
% self = basis({dof=0,description='',fftable=false, interval=[0,0]})
%
%     The generic basis class for all spectral expansions. This is a parent
%     class used by other classes. It is not meant to be useful as a
%     stand-alone class.
%
% Basis Properties:
%   domain - User-end domain after affine map from standard_domain
%   user_indexing - IndexingRule-derived object specifying a map from user-end indexing to the natural numbers 1, 2, ...
%   normalization - Normalization-derived object specifying the function normalizations
  properties
    description = [];
    fftable = false;
    normalization = []; % Normalization-derived object specifying the function normalizations
    user_indexing % IndexingRule-derived object specifying a map from user-end indexing to the natural numbers 1, 2, ...
    domain % { Interval1D([-1 1]) } | Other_Interval1D_object - Post-mapped univariate interval
  end
  properties(Abstract=true)
    standard_domain
    map_to_domain
    map_to_standard_domain
  end
  properties(SetAccess=protected)
    internal_indexing 
    internal_indexset
    allowed_function_normalizations
    default_function_normalization
    default_indexing_rule
  end
  properties(Access=protected)
    indexing_rules = IndexingRuleList.instance();
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
      %self.internal_indexing = self.indexing_parser(inp);
      self.internal_indexing = self.indexing_rules.parser(inp, self.default_indexing_rule);
      self.internal_indexset = self.internal_indexing.image;
    end
    function self = set.user_indexing(self,newindexing)
      %self.user_indexing = self.indexing_parser(newindexing);
      self.user_indexing = self.indexing_rules.parser(newindexing, self.default_indexing_rule);
    end
    function[self] = set.normalization(self, inp)
      self.normalization = self.function_normalization_parser(inp);
    end

    function output = range(self,N)
    % range -- Returns indices for the first N basis elements
    %
    % output = range(self,N)
    %
    %     Given an integer N, this function returns (user-end) indices for the
    %     first N basis elements.
      output = self.user_indexing(1:N);
    end

    function self = set.domain(self,newdomain)
      self.domain = newdomain;
      self.map_to_standard_domain = self.domain.compute_affine_map(self.standard_domain);
      self.map_to_domain = inv(self.map_to_standard_domain);
    end

    [n_array, nsize, numeln] = indexing(self,n);

  end
  methods(Access=protected)
    obj = function_normalization_parser(self, inp);
    %obj = indexing_parser(self, inp);
  end
end
