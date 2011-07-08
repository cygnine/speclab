classdef ZeroBasedIndexing < IndexingRule
  % ZeroBasedIndexing -- Zero-based indexing
  %
  % self = ZeroBasedIndexing()
  %
  %     Indexing that is zero-based: 0, 1, 2, ...
  properties(SetAccess=private)
    descriptive_adjective
    ids = {'0', 'zero', 'zeroes', 'zeros'};
    image
    dim = 1;
  end
  methods(Static)
    function self = instance()
      persistent obj
      if isempty(obj)
        obj = ZeroBasedIndexing();
      end
      self = obj;
    end
  end
  methods(Access=private)
    function self = ZeroBasedIndexing()
      self.descriptive_adjective = 'Zero-based';
      self.image = WholeNumbers.instance();
    end
  end
  methods
    function output = to_naturals(self, inp)
      output = inp + 1;
    end
    function output = from_naturals(self, inp)
      output = inp - 1;
    end
  end
end
