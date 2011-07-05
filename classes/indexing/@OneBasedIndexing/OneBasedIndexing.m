classdef OneBasedIndexing < IndexingRule
  % OneBasedIndexing -- One-based indexing
  %
  % self = OneBasedIndexing()
  %
  %     Indexing that is one-based: 1, 2, 3, ...
  properties(SetAccess=private)
    descriptive_adjective
    ids = {'1', 'one', 'ones'};
    image
  end
  methods(Static)
    function self = instance()
      persistent obj
      if isempty(obj)
        obj = OneBasedIndexing();
      end
      self = obj;
    end
  end
  methods(Access=private)
    function self = OneBasedIndexing()
      self.descriptive_adjective = 'One-based';
      self.image = NaturalNumbers.instance();
    end
  end
  methods
    function output = to_naturals(self, inp)
      output = inp;
    end
    function output = from_naturals(self, inp)
      output = inp;
    end
  end
end
