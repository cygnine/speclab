classdef IntegerBasis < Basis
% Integer -- Basis with "integer" indexing: 0, -1, 1, -2, 2, -3, 3, ...
%
% IntegerBasis(varargin)
% 
%     The basic class for spectral methods indexed 0, -1, 1, -2, 3, ..., varargin is
%     the same set of optional inputs as in the class Basis.
%
%     The purpose of this derived class is only for identification of the type
%     of indexing used for function identification in the basis -- derived basis
%     methods in general should not use the indexing properties of this class.
%     Rather, these indexing routines are meant to be used by external functions
%     to determine how the basis is indexed.
  methods
    function self = IntegerBasis(varargin)

      self = self@Basis(varargin{:});
    end

    function[ns] = range(N)
      persistent integer_range
      if isempty(integer_range)
        from speclab.common import integer_range
      end
      ns = integer_range(N);
    end

    subs = index_linear_to_natural(self, subs);
    subs = index_natural_to_linear(self, subs);
  end

end
