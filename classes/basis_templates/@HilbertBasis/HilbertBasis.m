classdef HilbertBasis < Basis
% HilbertBasis -- a class for orthogonal bases of Hilbert spaces
%
% self = HilbertBasis(...)
%
%     A generic (abstract) class derived from "Basis" classes -- HilbertBasis
%     instances correspond to complete bases of a Hilbert space whose elements 
%     are orthogonal with respect to the inner product defined on the Hilbert
%     space.
  properties(Access=protected)
    allowed_weight_normalizations = {};
    default_weight_normalization
  end
  properties
    weight_normalization = [];
  end
  methods
    function self = HilbertBasis(varargin)
      self = self@Basis(varargin{:});
    end
    function self = set.weight_normalization(self, inp)
      self.weight_normalization = self.weight_normalization_parser(inp);
    end
    obj = weight_normalization_parser(self, inp);
  end
  methods(Abstract=true)
    h = norm(self,n) % User-end method to compute Hilbert norms of basis element(s) n
  end
end
