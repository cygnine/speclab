classdef TotalGLexIndexing < IndexingRule
  % TotalGLexIndexing -- Total-degree-wise graded lexicographic indexing on R^N
  %
  % self = TotalGLexIndexing(dim)
  %
  %     The total ordering of N_0^dim corresponding to l^1 ('total') ordering, with
  %     intra-norm-value ordering done via lexicographic ordering.
  %
  %     The set N is the natural numbers.

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
      self.image = NaturalMultiIndices.instance(dim);
      self.dim = dim;
    end
  end
end
