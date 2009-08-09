classdef WholeBasis < basis.Basis
  methods
    function self = WholeBasis(varargin)
    % SELF = WHOLEBASIS(VARARGIN)
    % 
    %     The basic class for spectral methods indexed 0,1,2,.... VARARGIN is
    %     the same set of optional inputs as in the class Basis.
      self = basis.Basis(varargin{:});
    end
    function[ns] = range(N)
      ns = (0:(N-1)).';
    end
end
