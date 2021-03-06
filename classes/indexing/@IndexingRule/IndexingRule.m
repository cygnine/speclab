classdef IndexingRule < Singleton
% IndexingRule -- Superclass for indexing of bases
%
% self = IndexingRule()
%
%     A superclass for all singleton classes that are indexing rules. In
%     general, an indexing rule is a map from the natural numbers (1,2, ...) to
%     some prescribed index set. This map is (generally) an isomorphism. In
%     other words, this class endows sets corresponding to subclasses of
%     IndexSet with a total ordering.
%
%     All subclasses must define methods 'to_natural' and 'from_natural' that
%     implement the isomorphism to the natural numbers.
%
%     Generally speaking, the `forward' map is from the naturals to whatever the
%     custom range is. The 'backward' or 'inverse' maps takes the customized range
%     back to the naturals.
%
% IndexingRule Properties:
%   descriptive_adjective - Human-readable description of the indexing
%   ids - string and/or scalar identification for the indexing rule
%   image - An IndexSet that this rule maps to
%
% IndexingRule Methods:
%   range - Returns the first few indices
%   reindex - Maps indices from this rule into indices from another rule
%   to_naturals - Performs the map: from naturals to indices
%   from_naturals - Performs the inverse map: from indices to naturals
%   inv - Performs the inverse map: from indices to naturals
%   id_compare - Uses input id to identify if it identifies this class
  properties(SetAccess=private,Abstract=true)
    % e.g. 'one-based'
    descriptive_adjective
    ids % A cell array of strings that identify the rule.
    image % The IndexSet that is the range of this map
  end
  methods(Abstract=true)
    to_naturals
    from_naturals
  end
  properties(SetAccess=private)
    preimage = NaturalNumbers.instance(); % Always the naturals
  end
  methods
    function disp(self)
      fprintf('%s indexing\n', self.descriptive_adjective);
      inds = 1:10;
      for q = 1:length(inds)
        fprintf('[%2s] <-------> ', num2str(inds(q)));
        fprintf('[%s] \n', num2str(self.from_naturals(inds(q))));
      end
      fprintf('...\n');
    end
    function output = range(self, N)
    % range -- Returns the first N indices of the rule
    %
    % inds = range(self,N)
    % 
    %     When N is a natural number, this returns the first N indices of the
    %     indexing rule.
      output = self.from_naturals((1:N).');
    end
    function output = reindex(self, inp, other)
      % reindex -- Translates indices to another indexing
      %
      % output = reindex(self, inp, other)
      %
      %     Translates the input indices 'inp' that are from the indexing rule
      %     self into the corresponding indices of the indexing rule other.
      %     Since all IndexingRule's are isomorphic to the natural numbers, this
      %     is well-defined.

      if self==other
        output = inp;
      else
        output = other.from_naturals(self.to_naturals(inp));
      end
    end
    function output = inv(self,inp)
      % inv -- Inverse or `backward' map to naturals
      %
      % output = inv(self, inp)
      %
      %     Syntactic sugar for output = self.to_naturals(inp).
      output = self.to_naturals(inp);
    end
    function [tf, instance] = id_compare(self, inp)
    % id_compare -- Determines whether or not an input string ids this rule
    %
    % [tf, instance] = id_compare(self, inp)
    %
    %     Compares the string input 'inp' to all id strings in self.ids and
    %     returns true or false based on comparison.
      tf = any(strcmpi(num2str(inp), self.ids));
      if tf
        instance = self;
      else
        instance = [];
      end
    end
    function N = total_polynomial_space_dimension(self,d,k)
    % total_polynomial_space_dimension -- Returns the linear dimension of a polynomial space
    %
    % N = total_polynomial_space_dimension(d,k)
    %
    %     Returns the dimension of the space of d-variate polynomials of total
    %     degree not more than k. Is vectorized in the input k.

      N = zeros(size(k));

      for q = 1:numel(k);
        if k(q)<0
          N(q) = 0;
        else
          N(q) = nchoosek(d+k(q),d);
        end
      end
    end

  end
end
