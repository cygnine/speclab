function[x,w] = gauss_radau_quadrature(self,N,varargin)
% [x,w] = gauss_radau_quadrature(self,N, {r=self.domain.slices{1}.interval(2)})
%
%     Computes the n-point Gauss-Radau quadrature rule associated with the
%     system of orthogonal polynomials. If no optional input 'r' specifying the
%     location of the fixed point is given, it is assumed to be the right-hand
%     endpoint of the domain.

persistent strict_inputs gq
if isempty(strict_inputs)
  from labtools import strict_inputs
  from speclab.d1_utils import gauss_quadrature as gq
end

opt = strict_inputs({'r'}, {self.domain.slices{1}.interval(2)}, [], varargin{:});

[a,b] = self.recurrence(0:(N-1));
a = a(:); b = b(:);

temp = self.evaluate(opt.r, [N-2 N-1], 'normalization', MonicNormalization.instance());
%temp = eval_polynomial(r,a,b,[N-2,N-1], 'normalization','monic');
a(N) = opt.r - b(N)*temp(1)/temp(2);

%[x,w] = self.gauss_quadrature([], 'a', a, 'b', b);
[x,w] = gq(a,b);

x = self.map_to_domain(x);
w = self.scale_weight(w);
