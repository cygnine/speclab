classdef TensorBasis
% self = TensorBasis(Basis1, Basis2, ...)
%
%     Constructor for tensor-product basis. This results in a basis which has
%     ***** indexing over each dimension. Given N  Basis-type elements, this
%     constructs the tensor-product basis. This class assumes that the derived
%     Basis-types all have the following methods:
%
%      - .evaluate
%      - .quadrature_rule

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
