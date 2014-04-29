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
  % TensorProductBasis -- Tensor-product basis formed from univariate bases
  %
  % self = TensorProductBasis({Basis1, Basis2, Basis3, ...}, { key-value optional_inputs })
  %
  %     Given a cell array of univariate bases {Basis1, Basis2, ...}, this
  %     creates a basis with a tensor-product structure. The optional inputs
  %     indicate e.g. the structure of the user-end indexing.

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

    self = self@Basis(varargin{2:end});

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
    else
      %%%%%%%%%%%%%%%%% For special cases that needs fixing: 
      if isa(self.user_indexing, 'DegreeIndexing');
        self.user_indexing = DegreeIndexing.instance(self.dim);
      end
    end

    self.standard_domain = Orthotope('boundaries', standard_domains);
    self.domain = Orthotope('boundaries', domains);
    %self.map_to_standard_domain = self.domain.compute_affine_map(self.standard_domain);
    %self.map_to_domain = inv(self.map_to_standard_domain);
  end

  V = evaluate(self,x,n,varargin)
  w = weight(self, x, d);
  lmbda = derivative_expansion(self, n, d);
end

end
