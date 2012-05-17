classdef DegreeIndexing < IndexingRule
% DegreeIndexing -- Total-norm-wise degree indexing on N_0^dim
%
% self = DegreeIndexing(dim)
%
%     The total ordering of N_0^dim corresponding to l^1-norm ('total') ordering
%     -- thus this maps the naturals to the quotient space of N_0^dim formed by
%     defining all elements with an equivalent l^1 norm as an equivalence
%     class.
%
%     If one considers the non-quotient space N_0^dim, then this map is not an
%     isomorphism.
%
% DegreeIndexing Properties:
%   descriptive_adjective - Human-readable description of the indexing (here, 'l^1-norm-wise (''total'') degree indexing'')
%   ids - string and/or scalar identifications for this rule: 'degree'
%   image - The set WholeMultiIndices of the appropriate dimension
%
% DegreeIndexing Methods:
%   range - Returns the first few indices
%   reindex - Maps indices from this rule into indices from another rule
%   to_naturals - Performs the map: from naturals to indices
%   from_naturals - Performs the inverse map: from indices to naturals
%   inv - Performs the inverse map: from indices to naturals
%   id_compare - Returns this class if given 'degree'

  properties(SetAccess=private)
    descriptive_adjective
    ids = {'degree'};
    dim
    image
  end
  methods(Static)
    function self = instance(dim)
      persistent selfs
      if nargin==0
        dim = 1;
      else
      end
      if length(selfs)<dim
        self = DegreeIndexing(dim);
        selfs{dim} = self;
      else
        if isempty(selfs{dim})
          self = DegreeIndexing(dim);
          selfs{dim} = self;
        else
          self = selfs{dim};
        end
      end
    end
  end
  methods(Access=private)
    function self = DegreeIndexing(dim)
      self.descriptive_adjective = 'l^1-norm-wise (''total'') degree indexing';
      self.image = WholeMultiIndices.instance(dim);
      self.dim = dim;
    end
  end
end
