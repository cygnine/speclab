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
  properties(SetAccess=private)
    descriptive_adjective
    ids = {'integer'};
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
    function self = IntegerIndexing()
      self.descriptive_adjective = 'Integer';
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
