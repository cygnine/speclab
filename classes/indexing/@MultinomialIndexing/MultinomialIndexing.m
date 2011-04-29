classdef MultinomialIndexing < IndexingRule
  % MultinomialIndexing -- Indexing of the multinomials
  %
  % self = MultinomialIndexing(dim)
  %
  %     Indexing of the multinomial of dimension dim that proceeds:
  %
  %      0, -1, 1, -2, 2, -3, 3, -4, 4, ...
  %
  %     Most useful in cases like Fourier expansions.
  properties(SetAccess=private)
    descriptive_adjective
    dim
  end
  methods(Static)
    function self = instance(dim)
      persistent obj
      if length(obj)<dim
        obj = IntegerIndexing(dim);
      if isempty(obj{dim})
        self = IntegerIndexing(dim);
      else
        obj = selfs{dim};
      end
    end
  end
  methods(Access=Private)
    function self = MultinomialIndexing(dim)
      self.descriptive_adjective = 'Multinomial';
      if nargin==0
        self.dim = 1;
      end
    end
  end
  methods
    function output = to_naturals(self, inp)
      persistent linear_to_array_indexing
      if isempty(linear_to_array_indexing)
        from speclab.common.tensor import linear_to_array_indexing
      end
      error('Not implemented')
    end
    function output = from_naturals(self, inp)
      persistent linear_to_array_indexing
      if isempty(linear_to_array_indexing)
        from speclab.common.tensor import linear_to_array_indexing
      end
      output = linear_to_array_indexing(inp+1, dim);
    end
  end
end
