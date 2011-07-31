classdef MarginalGLexIndexing < IndexingRule
% MarginalGLexIndexing -- Total-degree-wise graded lexicographic indexing on R^N
%
% self = MarginalGLexIndexing(dim)
%
%     The total ordering of N_0^dim corresponding to marginal-degree ordering, with
%     intra-degree ordering done via (standard) lexicographic ordering.
%
%     The set N_0 is the naturals with 0: {0, 1, 2, ... }
%
% MarginalGLexIndexing Properties:
%   descriptive_adjective - Human-readable description of the indexing (here, 'Marginal-dimension-wise graded lexicographic')
%   ids - string and/or scalar identifications for this rule: 'marginalglex', 'mglex'
%   image - The set WholeMultiIndices of the appropriate dimension
%
% MarginalGLexIndexing Methods:
%   range - Returns the first few indices
%   reindex - Maps indices from this rule into indices from another rule
%   to_naturals - Performs the map: from naturals to indices
%   from_naturals - Performs the inverse map: from indices to naturals
%   inv - Performs the inverse map: from indices to naturals
%   id_compare - Returns this class if given 'marginalglex' or 'mglex'
  properties(SetAccess=private)
    descriptive_adjective
    ids = {'marginalglex', 'mglex'};
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
        self = MarginalGLexIndexing(dim);
        selfs{dim} = self;
      else
        if isempty(selfs{dim})
          self = MarginalGLexIndexing(dim);
          selfs{dim} = self;
        else
          self = selfs{dim};
        end
      end
    end
  end
  methods(Access=private)
    function self = MarginalGLexIndexing(dim)
      self.descriptive_adjective = 'Marginal-dimension-wise graded lexicographic';
      self.image = WholeMultiIndices.instance(dim);
      self.dim = dim;
    end
  end
end
