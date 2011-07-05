classdef TotalGRevLexIndexing < IndexingRule
% TotalGRevLexIndexing -- Total-degree-wise reverse lexicographic indexing on R^N
%
% self = TotalGRevLexIndexing(dim)
%
%     The total ordering of N_0^dim corresponding to total-degree ordering, with
%     intra-degree ordering done via reverse lexicographic ordering.
%
%     The set N_0 is the naturals with 0: {0, 1, 2, ... }

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
          self = TotalGRevlexIndexing(dim);
          selfs{dim} = self;
        else
          self = selfs{dim};
        end
      end
    end
  end
  methods(Access=private)
    function self = TotalGRevLexIndexing(dim)
      self.descriptive_adjective = 'Total-degree-wise graded reverse lexicographic';
      self.image = WholeMultiIndices.instance(dim);
      self.dim = dim;
    end
  end
end
