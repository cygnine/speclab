function[x,w] = gauss_lobatto_quadrature(self,N,varargin)
% [x,w] = gauss_radau_quadrature(self,N, {r1=self.domain.slices{1}.interval(1), 
%                                         r2=self.domain.slices{1}.interval(2)})
%
%     Computes the n-point Gauss-Lobatto quadrature rule associated with the
%     system of orthogonal polynomials. If no optional inputs 'r1' and 'r2'
%     specifying the location of the fixed points are given, they are assumed to
%     be the endpoints of the domain.

persistent parser input_parser
if isempty(parser)
  from labtools import input_parser
  [opt, parser] = input_parser({'r1', 'r2'}, ...
                   {self.standard_domain.interval(1), self.standard_domain.interval(2)}, ...
                   [], ...
                   varargin{:});
else
  parser.parse(varargin{:});
  opt = parser.Results;
end

% Call standard Gauss quadrature driver with modified terminal recurrence
% coefficient.
[a,b] = self.recurrence(0:(N-1));
a = a(1:N);
b = b(1:N);

domain_storage = self.domain;
self.domain = self.standard_domain;
inds = self.range(N); inds = inds([end end-1]);
temp = self.evaluate([opt.r1; opt.r2], inds, 'normalization', MonicNormalization.instance());
self.domain = domain_storage;

modif = inv(temp)*[opt.r1*temp(1,1); opt.r2*temp(2,1)];

% Lobatto modification for Jacobi matrix: 1999_gautschi
a(N) = modif(1);
b(N) = modif(2);

[x,w] = OrthogonalPolynomialBasis.gauss_quadrature_driver(a,b);
[x,w] = self.scale_quadrature(x,w);
