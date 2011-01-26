function[x,w] = gauss_lobatto_quadrature(self,N,varargin)
% [x,w] = gauss_radau_quadrature(self,N, {r1=self.domain.slices{1}.interval(1), 
%                                         r2=self.domain.slices{1}.interval(2)})
%
%     Computes the n-point Gauss-Lobatto quadrature rule associated with the
%     system of orthogonal polynomials. If no optional inputs 'r1' and 'r2'
%     specifying the location of the fixed points are given, they are assumed to
%     be the endpoints of the domain.

persistent strict_inputs gq
if isempty(gq)
  from labtools import strict_inputs
  from speclab.d1_utils import gauss_quadrature as gq
end
opt = strict_inputs({'r1', 'r2'}, {self.domain.slices{1}.interval(1), ...
                                   self.domain.slices{1}.interval(2)}, [], varargin{:});

% Move all inputs to the standard interval
opt.r1 = self.map_to_standard_domain(opt.r1);
opt.r2 = self.map_to_standard_domain(opt.r2);

[a,b] = self.recurrence(0:(N-1));
a = a(1:N);
b = b(1:N);

temp = self.evaluate([opt.r1; opt.r2],[N-1,N-2],'normalization',MonicNormalization.instance());
modif = inv(temp)*[opt.r1*temp(1,1); opt.r2*temp(2,1)];

% Lobatto modification for Jacobi matrix: 1999_gautschi
a(N) = modif(1);
b(N) = modif(2);

[x,w] = gq(a,b);

% Convert back to physical inteval
x = self.map_to_domain(x);
w = self.scale_weight(w);
