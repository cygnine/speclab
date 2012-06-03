function[x,w] = gauss_radau_quadrature(self,N,varargin)
% [x,w] = gauss_radau_quadrature(self,N, {r=self.domain.interval(2)})
%
%     Computes the n-point Gauss-Radau quadrature rule associated with the
%     system of orthogonal polynomials. If no optional input 'r' specifying the
%     location of the fixed point is given, it is assumed to be the right-hand
%     endpoint of the domain.

persistent parser input_parser
if isempty(parser)
  from labtools import input_parser
  [opt, parser] = input_parser({'r'}, ...
                        {self.map_to_standard_domain(self.domain.interval(2))}, ...
                        [],...
                        varargin{:});
else
  parser.parse(varargin{:});
  opt = parser.Results;
end

[a,b] = self.recurrence(0:(N-1));

domain_storage = self.domain;
self.domain = self.standard_domain;
inds = self.range(N); inds = inds([end-1 end]);
%temp = self.evaluate(opt.r, inds, 'normalization', 'monic');
temp = self.evaluate(opt.r, inds, 'normalization', MonicNormalization.instance());
a(N) = opt.r - b(N)*temp(1)/temp(2);
self.domain = domain_storage;

[x,w] = OrthogonalPolynomialBasis.gauss_quadrature_driver(a,b);
[x,w] = self.scale_quadrature(x,w);
