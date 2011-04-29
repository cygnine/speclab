classdef TensorProductBasis
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
    
    self.bases = {};
    for q = 1:nargin
      if isa(varargin{q}, 'Basis');
        self.bases{end+1} = varargin{q};
      end
    end
    self.dim = length(self.bases);

    self.indexing = indexing;
  end

  V = evaluate(self,x,n)
end

end
