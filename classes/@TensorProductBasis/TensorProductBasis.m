classdef TensorProductBasis < Basis
properties(SetAccess=protected)
  bases
  dim
  indexing
end
methods 
  function[self] = TensorProductBasis(varargin)

    persistent indexing
    if isempty(indexing)
      from speclab.common.tensor import linear_to_array_indexing as indexing
    end

    self = self@Basis(varargin{:});

    self.bases = {};
    if nargin == 2
      % If user said "TensorProductBasis(SomeBasis, someinteger)"
      if isa(varargin{2}, 'numeric')
        self.bases{1} = varargin{1};
        self.bases{2} = varargin{1};
      end
    end
    % If user said "TensorProductBasis(SomeBasis, OtherBasis, ...)"
    if isempty(self.bases);
      for q = 1:nargin
        if isa(varargin{q}, 'Basis');
          self.bases{end+1} = varargin{q};
        end
      end
    end
    self.dim = length(self.bases);

    self.indexing = indexing;
  end

  V = evaluate(self,x,n)
end

end
