classdef MarginalGRevLexIndexing < IndexingRule
% MarginalGRevLexIndexing -- Marginal-dimension-wise reverse lexicographic indexing on R^N
%
% self = MarginalGRevLexIndexing(dim)
%
%     The marginal ordering of N^dim corresponding to l^\inf norm ('marginal')
%     ordering, with intra-norm-value ordering done via reverse lexicographic
%     ordering.
%
%     The set N is the natural numbers.
%
% MarginalGRLexIndexing Properties:
%   descriptive_adjective - Human-readable description of the indexing (here, 'Marginal-dimension-wise graded reverse lexicographic')
%   ids - string and/or scalar identifications for this rule: 'marginalgrevlex', 'mgrevlex'
%   image - The set WholeMultiIndices of the appropriate dimension
%
% MarginalGRLexIndexing Methods:
%   range - Returns the first few indices
%   reindex - Maps indices from this rule into indices from another rule
%   to_naturals - Performs the map: from naturals to indices
%   from_naturals - Performs the inverse map: from indices to naturals
%   inv - Performs the inverse map: from indices to naturals
%   id_compare - Returns this class if given 'marginalgrevlex' or 'mgrevlex'
  properties(SetAccess=private)
    descriptive_adjective
    ids = {'marginalgrevlex', 'mgrevlex'};
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
        self = MarginalGRevLexIndexing(dim);
        selfs{dim} = self;
      else
        if isempty(selfs{dim})
          self = MarginalGRevlexIndexing(dim);
          selfs{dim} = self;
        else
          self = selfs{dim};
        end
      end
    end
  end
  methods(Access=private)
    function self = MarginalGRevLexIndexing(dim)
      self.descriptive_adjective = 'Marginal-dimension-wise graded reverse lexicographic';
      self.image = WholeMultiIndices.instance(dim);
      self.dim = dim;
    end
  end
end
