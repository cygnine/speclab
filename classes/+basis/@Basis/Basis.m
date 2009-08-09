classdef Basis
  properties
    dof = 0;
    N = 0;
    description = 'None';
    fftable = false;
    vandermonde = [];
    vandermonde_inverse = [];
    evaluate = [];
    derivative = [];
    range = [];
    differentiation_matrix = [];
    interval = [0,0];
    work_functions = false;
    quadrature_rule = [];
  end
  methods
    function self = Basis(varargin)
    % SELF = BASIS({DOF=0,DESCRIPTION='None',FFTABLE=false})
    %
    %     The generic basis class for all spectral expansions.

      global handles;
      inputs = {'dof','description','fftable'};
      defaults = {0, 'None',false};
      self = handles.common.InputSchema(inputs, defaults, [], varargin{:});
    end
    
    function value = get.N(self)
      value = self.dof;
    end
  end
end
