classdef JacobiQuadratureRule < QuadratureRule
  properties
  end
  methods
    function self = JacobiQuadratureRule(dof,varargin);
      % self = JacobiQuadratureRule(dof, {type='gauss',alpha=-1/2, beta=-1/2, scale=1, shift=0, r=1})
      %
      %     Creates a JacobiQuadratureRule object. dof is the total number of
      %     degrees of freedom. The optional inputs are the usual specifiers for
      %     Jacobi Polynomial family objects. 

      if dof<1
        error('The degrees of freedom must be a positive integer');
      end
      global handles;
      jac = handles.speclab.orthopoly1d.jacobi;
      inputs = {'type', 'alpha', 'beta', 'interval', 'scale', 'shift', 'r'};
      defaults = {'gauss', -1/2, -1/2, [], 1, 0, 1};
      opt = handles.common.InputSchema(inputs, defaults, [], varargin{:});
      [params.alpha, params.beta, params.scale, params.shift] = deal(...
        opt.alpha, opt.beta, opt.scale, opt.shift);


      if strcmpi(opt.type, 'gauss')
        [nodes, weights] = jac.quad.gauss_quadrature(dof,params);
        descr = 'Jacobi-Gauss quadrature rule';
      elseif strcmpi(opt.type, 'gauss-radau')
        [nodes, weights] = ...
          jac.quad.gauss_radau_quadrature(dof, params);
        descr = 'Jacobi-Gauss-Radau quadrature rule';
      elseif strcmpi(opt.type, 'gauss-lobatto')
        [nodes, weights] = ...
          jac.quad.gauss_lobatto_quadrature(dof, params);       
        descr = 'Jacobi-Gauss-Lobatto quadrature rule';
      else
        error('I don''t recognize the input argument ''type''');
      end
      self = self@QuadratureRule('dof', dof, 'nodes', nodes, ...
                                 'weights', weights, 'description', descr);
      self.parameters = params;
    end

    function vals = weight_function(self,x)
      global handles;
      jac = handles.speclab.orthopoly1d.jacobi;
      vals = jac.weights.weight(x,self.parameters);
    end
end
