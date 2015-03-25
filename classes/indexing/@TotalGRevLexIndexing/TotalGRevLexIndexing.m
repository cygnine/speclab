classdef TotalGRevLexIndexing < IndexingRule
% TotalGRevLexIndexing -- Total-norm-wise reverse lexicographic indexing on N^dim
%
% self = TotalGRevLexIndexing(dim)
%
%     The total ordering of N^dim corresponding to l^1-norm ('total') ordering, with
%     intra-norm-value ordering done via reverse lexicographic ordering.
%
%     The set N is the natural numbers.
%
% TotalGRevLexIndexing Properties:
%   descriptive_adjective - Human-readable description of the indexing (here, 'l^1-norm-wise (''total'') graded reverse lexicographic'')
%   ids - string and/or scalar identifications for this rule: 'totalgrevlex', 'tgrevlex'
%   image - The set WholeMultiIndices of the appropriate dimension
%
% TotalGRevLexIndexing Methods:
%   range - Returns the first few indices
%   reindex - Maps indices from this rule into indices from another rule
%   to_naturals - Performs the map: from naturals to indices
%   from_naturals - Performs the inverse map: from indices to naturals
%   inv - Performs the inverse map: from indices to naturals
%   id_compare - Returns this class if given 'totalgrevlex' or 'tgrevlex'

  properties(SetAccess=private)
    descriptive_adjective
    ids = {'totalgrevlex', 'tgrevlex'};
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
        self = TotalGRevLexIndexing(dim);
        selfs{dim} = self;
      else
        if isempty(selfs{dim})
          self = TotalGRevLexIndexing(dim);
          selfs{dim} = self;
        else
          self = selfs{dim};
        end
      end
    end
  end
  methods(Access=private)
    function self = TotalGRevLexIndexing(dim)
      self.descriptive_adjective = 'l^1-norm-wise (''total'') graded reverse lexicographic';
      self.image = WholeMultiIndices.instance(dim);
      self.dim = dim;
    end
  end
end
