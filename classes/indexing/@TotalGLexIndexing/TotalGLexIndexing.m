classdef TotalGLexIndexing < IndexingRule
% TotalGLexIndexing -- Total-degree-wise graded lexicographic indexing on R^N
%
% self = TotalGLexIndexing(dim)
%
%     The total ordering of N_0^dim corresponding to l^1 ('total') ordering, with
%     intra-norm-value ordering done via lexicographic ordering.
%
%     The set N is the natural numbers.
%
% TotalGLexIndexing Properties:
%   descriptive_adjective - Human-readable description of the indexing (here, 'l^1-norm-wise (''total'') graded lexicographic'')
%   ids - string and/or scalar identifications for this rule: 'totalglex', 'tglex'
%   image - The set WholeMultiIndices of the appropriate dimension
%
% TotalGLexIndexing Methods:
%   range - Returns the first few indices
%   reindex - Maps indices from this rule into indices from another rule
%   to_naturals - Performs the map: from naturals to indices
%   from_naturals - Performs the inverse map: from indices to naturals
%   inv - Performs the inverse map: from indices to naturals
%   id_compare - Returns this class if given 'totalglex' or 'tglex'

  properties(SetAccess=private)
    descriptive_adjective
    ids = {'totalglex', 'tglex'};
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
        self = TotalGLexIndexing(dim);
        selfs{dim} = self;
      else
        if isempty(selfs{dim})
          self = TotalGLexIndexing(dim);
          selfs{dim} = self;
        else
          self = selfs{dim};
        end
      end
    end
  end
  methods(Access=private)
    function self = TotalGLexIndexing(dim)
      self.descriptive_adjective = 'l^1-norm-wise (''total'') graded lexicographic';
      self.image = WholeMultiIndices.instance(dim);
      self.dim = dim;
    end
  end
end
