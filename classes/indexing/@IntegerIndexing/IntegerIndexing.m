classdef IntegerIndexing < IndexingRule
% IntegerIndexing -- Indexing of the integers
%
% self = IntegerIndexing()
%
%     Indexing of the integer that proceeds:
%
%      0, -1, 1, -2, 2, -3, 3, -4, 4, ...
%
%     Most useful in cases like Fourier expansions.
%
% IntegerIndexing Properties:
%   descriptive_adjective - Human-readable description of the indexing (here, 'Integer')
%   ids - string expressions identifying this rule: 'int', 'integer'
%   image - The set IntegerNumbers
%
% IntegerIndexing Methods:
%   range - Returns the first few indices: 0, -1, 1, ...
%   reindex - Maps indices from this rule into indices from another rule
%   to_naturals - Performs the map: from naturals to integers
%   from_naturals - Performs the inverse map: from integers to naturals
%   inv - Performs the inverse map: from integers to naturals
%   id_compare - Returns this class if given 'integer' or 'int'
  properties(SetAccess=private)
    descriptive_adjective
    ids = {'integer', 'int'};
    image
  end
  methods(Static)
    function self = instance()
      persistent obj
      if isempty(obj)
        obj = IntegerIndexing();
      end
      self = obj;
    end
  end
  methods(Access=private)
    function self = IntegerIndexing(varargin)
      self.descriptive_adjective = 'Integer';
      self.image = IntegerNumbers.instance();
    end
  end
  methods
    function output = to_naturals(self, inp)
      signs = (sign(inp)-1)/2;
      output = 1 + 2*abs(inp) + signs;
      output(inp==0) = 1;
    end
    function output = from_naturals(self, inp)
      signs = mod(inp,2)*2 - 1;
      output = signs.*(ceil((inp-1)/2));
    end
  end
end
