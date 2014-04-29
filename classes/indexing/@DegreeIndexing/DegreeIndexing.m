classdef DegreeIndexing < IndexingRule
% DegreeIndexing -- Total-norm-wise degree indexing on N_0^dim
%
% self = DegreeIndexing(dim)
%
%     Maps the naturals to the set of whole numbers (naturals plus 0). The
%     degree of a dim-dimensional multi-index (on N_0^dim) is defined as its
%     l^1 norm. There is 1 multi-index with degree 0. There are dim
%     multi-indices with degree 1, etc. This mapping effects the connection
%     between a degree-based ordering of the basis and the polynomial degree.
%
%     Thus, for degree 2:
%
%        natural ordering         degree 
%            1       <--------->     0
%            2       <--------->     1
%            3       <--------->     1
%            4       <--------->     2
%            5       <--------->     2
%            6       <--------->     2
%
%     This map is not an isomorphism, except when dim = 1. When mapping *to*
%     naturals, this mapping is one-to-many. When mapping *from* naturals, it's
%     many-to-one.
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
    ids = {'deg', 'degree'};
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
      self.descriptive_adjective = 'degree indexing based on N_0^dim ordering';
      self.image = WholeNumbers.instance();
      self.dim = dim;
    end
  end
end
