classdef TotalGRevLexIndexing < IndexingRule
% TotalGRevLexIndexing -- Total-norm-wise reverse lexicographic indexing on R^N
%
% self = TotalGRevLexIndexing(dim)
%
%     The total ordering of N^dim corresponding to l^1-norm ('total') ordering, with
%     intra-norm-value ordering done via reverse lexicographic ordering.
%
%     The set N is the natural numbers.

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
      self.descriptive_adjective = 'l^1-norm-wise (''total'') graded reverse lexicographic';
      self.image = NaturalMultiIndices.instance(dim);
      self.dim = dim;
    end
  end
end
