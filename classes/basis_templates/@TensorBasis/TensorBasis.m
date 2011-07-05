classdef TensorBasis
% self = TensorBasis(Basis1, Basis2, ...)
%
%     Constructor for tensor-product basis. The indexing is given by some
%     variant of a multidimensional indexing (e.g. TotalGRevLexIndexing). Given
%     N  Basis-type elements, this constructs the tensor-product basis. This
%     class assumes that the derived Basis-types all have the following methods:
%
%      - .evaluate

  properties(Access=private)
    bases = {};
  end
  methods 
    function self = TensorBasis(varargin)
    
      N = nargin;
      for q = 1:nargin
        if not(isa(varargin{q}, 'Basis'));
          error('All inputs must be of Basis-type');
        end
      end

      self.bases = varargin;
    end
  end

end
