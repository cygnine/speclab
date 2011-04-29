classdef WholeBasis < Basis
  methods
    function self = WholeBasis(varargin)
    % WholeBasis -- Basis with 'natural' indexing based on whole numbers 0, 1, ...
    %
    % WholeBasis(varargin)
    % 
    %     The basic class for spectral methods indexed 0,1,2,.... VARARGIN is
    %     the same set of optional inputs as in the class Basis.
    %
    %     The purpose of this derived class is only for identification of the type
    %     of indexing used for function identification in the basis -- derived basis
    %     methods in general should  not use the indexing properties of this class.
    %     Rather, these indexing routines are meant to be used by external functions
    %     to determine how the basis is indexed.

      self = self@Basis(varargin{:});
    end

    function[ns] = range(N)
      ns = (0:(N-1)).';
    end

    subs = natural_indexing(self, subs);
    subs = index_natural_to_linear(self, subs);
  end
end
