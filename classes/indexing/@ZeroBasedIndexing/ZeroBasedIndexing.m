classdef ZeroBasedIndexing < IndexingRule
% ZeroBasedIndexing -- Zero-based indexing
%
% self = ZeroBasedIndexing()
%
%     Indexing that is zero-based: 0, 1, 2, ...
%
% ZeroBasedIndexing Properties:
%   descriptive_adjective - Human-readable description of the indexing (here, 'Zero-based')
%   ids - string expressions identifying this rule: 0, 'zero', 'zeroes', 'zeros'
%   image - The set WholeNumbers
%
% ZeroBasedIndexing Methods:
%   range - Returns the first few indices: 0, 1, 2, ...
%   reindex - Maps indices from this rule into indices from another rule
%   to_naturals - Performs the map: from naturals to whole numbers
%   from_naturals - Performs the inverse map: from whole numbers to naturals
%   inv - Performs the inverse map: from whole numbers to naturals
%   id_compare - Returns this class if given 0, 'zero', 'zeroes', or 'zeros'
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
