function[x,w] = gauss_quadrature(self,n)
% [x,w] = gauss_quadrature(self,n)
%
%     Computes the n-point Gauss quadrature rule associated with the system of
%     orthogonal polynomials. 

persistent gq
if isempty(gq)
  from speclab.d1_utils import gauss_quadrature as gq
end
%opt = strict_inputs({'a', 'b'}, {[], []}, [], varargin{:});

%if isempty(opt.a)
%  J = self.jacobi_matrix(n);
%else
%  opt.a = a(:); opt.b = b(:);
%  n = min([length(opt.a), length(opt.b)]);
%  J = spdiags([sqrt([b(2:n);0]) a(1:n) sqrt(b(1:n))], -1:1, n, n);
%end
%
%x = eig(J);
%w = 1./sum(self.evaluate(x, 0:(n-1), 'normalization', ...
%                          OrthonormalNormalization.instance()).^2,2);
%

[a,b] = self.recurrence(0:(n-1));
[x,w] = gq(a, b);

x = self.map_to_domain(x);
w = self.scale_weight(w);
