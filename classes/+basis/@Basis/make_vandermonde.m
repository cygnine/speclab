function self = make_vandermonde(self)
% SELF = MAKE_VANDERMONDE(SELF)
%
%     Uses the data stored in the class instance SELF to construct an
%     interpolating Vandermonde matrix. This method is not valid in this base
%     class. 
  
% First determine what nodal set for interpolation to use.
if isempty(self.quadrature_rule)
  if isempty(self.interpolation_nodes)
    self = self.set_default_quadrature_rule;
    self.vandermonde = self.evaluate(self.quadrature_rule.nodes, ...
       self.range(self.dof));
    self.vandermonde_inverse = self.vandermonde'*...
      spdiags(self.quadrature_rule.weights,0,self.dof,self.dof);
  elseif length(self.interpolation_nodes)~=self.dof
    error('You can''t have a different # of interpolation nodes and degrees of freedom unless you give me a quadrature rule');
  else
    self.vandermonde = self.evaluate(self.interpolation_nodes,...
      self.range(self.dof));
  end
else
  self.vandermonde = self.evaluate(self.quadrature_rule.nodes, ...
      self.range(self.dof));
end
