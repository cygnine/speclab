function[x,w] = gauss_quadrature(self,n)
% [x,w] = gauss_quadrature(self,n)
%
%     Computes the n-point Gauss quadrature rule associated with the system of
%     Chebyshev polynomials.

assert(n > 0);
assert(numel(n) == 1, 'Input n must be scalar');

k = (0:(n-1)).';

x = flipud(cos(pi*(2*k+1)/(2*n)));
w = pi/n*ones(size(x));

[x,w] = self.scale_quadrature(x,w);
