classdef QuadratureRule
  properties
    dof = 0;
    description = 'no description');
    nodes = [];
    weights = [];
    parameters = [];
  end
  methods
    function self = QuadratureRule(varargin);
      % self = QuadratureRule({dof=0, nodes = [], weights = [], 
      %                        description = 'no description'})
      %
      %     The base class instance for quadrature rules.

      global packages;
      inputs = {'dof', 'nodes', 'weights', 'description'};
      defaults = {0,[],[],'no description'};
      self = packages.labtools.input_schema(inputs,defaults,[],varargin{:});
    end
    function vals = weight_function(self,x)
      vals = 0*x;
    end
  end
end
