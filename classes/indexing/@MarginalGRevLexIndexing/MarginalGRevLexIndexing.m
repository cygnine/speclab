classdef MarginalGRevLexIndexing < IndexingRule
  % MarginalGRevLexIndexing -- Marginal-degree-wise reverse lexicographic indexing on R^N
  %
  % self = MarginalGRevLexIndexing(dim)
  %
  %     The total ordering of N_0^dim corresponding to total-degree ordering, with
  %     intra-degree ordering done via reverse lexicographic ordering.
  %
  %     The set N_0 is the naturals with 0: {0, 1, 2, ... }

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
      self.descriptive_adjective = 'Marginal-degree-wise graded reverse lexicographic';
      self.image = WholeMultiIndices.instance(dim);
      self.dim = dim;
    end
  end
end
