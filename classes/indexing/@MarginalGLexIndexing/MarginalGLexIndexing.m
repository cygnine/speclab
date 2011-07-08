classdef MarginalGLexIndexing < IndexingRule
  % MarginalGLexIndexing -- Total-degree-wise graded lexicographic indexing on R^N
  %
  % self = MarginalGLexIndexing(dim)
  %
  %     The total ordering of N_0^dim corresponding to marginal-degree ordering, with
  %     intra-degree ordering done via (standard) lexicographic ordering.
  %
  %     The set N_0 is the naturals with 0: {0, 1, 2, ... }

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
      self.image = NaturalMultiIndices.instance(dim);
      self.dim = dim;
    end
  end
end
