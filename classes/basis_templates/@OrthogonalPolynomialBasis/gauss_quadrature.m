function[x,w] = gauss_quadrature(self,n)
% [x,w] = gauss_quadrature(self,n)
%
%     Computes the n-point Gauss quadrature rule associated with the system of
%     orthogonal polynomials. 

[a,b] = self.recurrence(0:(n-1));
[x,w] = OrthogonalPolynomialBasis.gauss_quadrature_driver(a,b);
[x,w] = self.scale_quadrature(x,w);
