classdef DirectSumIndexingRule < IndexingRule
% DirectSumIndexingRule -- Glue to index the direct sum of one-dimensional IndexingRule's
%
% self = DirectSumIndexingRule(MultiDimIndexingRule, UnivarRule1, UnivarRul2, ...)
%
%     Let d be the number of UnivarRule's provided. This class provides a way to
%     enumerate the direct sum of the images of UnivarRule's. Concretely,
%     MultiDimIndexingRule must satisfy:
%       (a) MultiDimIndexingRule.dim == d 
%       (b) MultiDimIndexingRule.image is the NaturalMultiIndices of d dimensions
%
%     With this, then MultiDimIndexingRule maps N ---> N^d. Then for any j \in
%     N, let z = MultiDimIndexingRule(j). This rule will then return y, a vector
%     of length d, satisfying y(m) = UnivarRule(z(m)) for m = 1, 2, ..., d.
  properties(SetAccess=private)
    % e.g. 'one-based'
    descriptive_adjective
    ids % A cell array of strings that identify the rule.
    image % The IndexSet that is the range of this map
  end
  properties(SetAccess=private)
    rules = {};
    multirule;
    dim;
  end
  methods(Static)
    function[self] = instance(varargin);
    % This is a wrapper to return essentially an empty rule
      self = DirectSumIndexingRule();
    end
  end
  methods
    function[self] = DirectSumIndexingRule(multidimrule, varargin)
      self.descriptive_adjective = 'Direct sum';
      self.ids = {};
      self.image = [];

      if nargin == 0 
        % Return an empty rule
        return;
      end
      d = length(varargin);
      if d < 1
        error('At least two inputs must be given');
      elseif not(isa(multidimrule, 'IndexingRule'))
        error('First input must be an IndexingRule object');
      elseif (d ~= multidimrule.dim) & (d > 1)
        error('The dimension of the image of the first rule must be the number of remaining rules')
      end
      for q = 1:d
        if not(isa(varargin{q}, 'IndexingRule'))
          error('All inputs must be IndexingRule objects');
        end
      end

      % Down to business:
      self.dim = multidimrule.dim;
      self.multirule = multidimrule;
      if d == 1
        for q = 1:self.dim
          self.rules{q} = varargin{1};
        end
      else
        self.rules = varargin;
      end
    end
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
    n = to_naturals(self, a);
    a = from_naturals(self, n);
  end
end
