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
      self.image = NaturalMultiIndices.instance(dim);
      self.dim = dim;
    end
  end
end
