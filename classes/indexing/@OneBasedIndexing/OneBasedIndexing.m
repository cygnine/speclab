classdef OneBasedIndexing < IndexingRule
% OneBasedIndexing -- One-based indexing
%
% self = OneBasedIndexing()
%
%     Indexing that is one-based: 1, 2, 3, ...
%
% OneBasedIndexing Properties:
%   descriptive_adjective - Human-readable description of the indexing (here, 'One-based')
%   ids - string expressions identifying this rule: 1, 'one', 'ones'
%   image - The set NaturalNumbers
%
% OneBasedIndexing Methods:
%   range - Returns the first few indices: 1, 2, ...
%   reindex - Maps indices from this rule into indices from another rule
%   to_naturals - Performs the map: the identity
%   from_naturals - Performs the inverse map: the identity
%   inv - Performs the inverse map: the identity
%   id_compare - Returns this class if given 1, 'one', or 'ones'
  properties(SetAccess=private)
    descriptive_adjective
    ids = {'1', 'one', 'ones'};
    image
    dim = 1;
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
