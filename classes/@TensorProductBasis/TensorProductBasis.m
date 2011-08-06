classdef TensorProductBasis < Basis
properties(SetAccess=protected)
  bases
  dim
  %indexing
end
properties
  standard_domain
  map_to_domain
  map_to_standard_domain
end
methods 
  function[self] = TensorProductBasis(varargin)

    %persistent indexing
    %if isempty(indexing)
    %  from speclab.common.tensor import linear_to_array_indexing as indexing
    %end

    persistent inparse
    if isempty(inparse)
      inparse = inputParser();
      inparse.KeepUnmatched = true;

      inparse.addRequired('bases', @(x) isa(x, 'cell'));
      inparse.addParamValue('indexing', []);
      inparse.addParamValue('internal_indexing', []);
    end

    inparse.parse(varargin{:});
    parsed_inputs = inparse.Results;
    bases = parsed_inputs.bases;

    %inputs = {'indexing', 'internal_indexing'};
    %defaults = {[], []};

    self = self@Basis(varargin{:});

    self.bases = {};
    internal_rules = {};
    user_rules = {};
    standard_domains = {};
    domains = {};

    if nargin == 2
      % If user said "TensorProductBasis(SomeBasis, someinteger)"
      if isa(bases{2}, 'numeric')
        for q = 1:bases{2}
          self.bases{q} = bases{1};
          internal_rules{q} = self.bases{q}.internal_indexing;
          user_rules{q} = self.bases{q}.user_indexing;
          standard_domains{q} = self.bases{q}.standard_domain;
          domains{q} = self.bases{q}.domain;
        end
      end
    end
    % If user said "TensorProductBasis(SomeBasis, OtherBasis, ...)"
    if isempty(self.bases);
      for q = 1:length(bases)
        if isa(bases{q}, 'Basis');
          self.bases{q} = bases{q};
          internal_rules{q} = self.bases{q}.internal_indexing;
          user_rules{q} = self.bases{q}.user_indexing;
          standard_domains{q} = self.bases{q}.standard_domain;
          domains{q} = self.bases{q}.domain;
        end
      end
    end
    self.dim = length(self.bases);

    self.default_indexing_rule = TotalGRevLexIndexing.instance(self.dim);
    self.internal_indexing = parsed_inputs.internal_indexing;
    self.user_indexing = parsed_inputs.indexing;

    self.internal_indexing = DirectSumIndexingRule(self.internal_indexing, internal_rules{:});
    if isempty(parsed_inputs.indexing)
      self.user_indexing = DirectSumIndexingRule(self.user_indexing, user_rules{:});
    end

    self.standard_domain = Orthotope('boundaries', standard_domains);
    self.domain = Orthotope('boundaries', domains);
    %self.map_to_standard_domain = self.domain.compute_affine_map(self.standard_domain);
    %self.map_to_domain = inv(self.map_to_standard_domain);
  end

  V = evaluate(self,x,n)
end

end
