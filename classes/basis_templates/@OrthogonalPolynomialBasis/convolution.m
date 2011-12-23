function[w] = convolution(self, u, v, other)
% convolution -- Discrete convolution of spectral expansions 
%
% w = convolution(self, u, v, [[ other = self ]])
%
%     u is a length-N coefficient vector for a polynomial expansion. v is a
%     length-M coefficient vector. The output w is a length-(N+M-1) coefficient
%     vector corresponding to the product of the polynomials that u and v
%     represent.
%
%     u, and w represent coefficients corresponding to the normalization of
%     self. If the input "other" is not given, then the same is assumed for v.
%     If "other" is a given OrthogonalPolynomialBasis instance, v is a
%     coefficient vector for a polynomial in that expansion.

if nargin > 3
  error('Not yet implemented for other inputs'); 
end

N = length(u);
M = length(v);

w = zeros([N + M - 1 1]);

% Need triple product arrays:
ek = self.triple_product(max([N M]), N+M-1);

% This computes inner products of the functional product
for q = 1:(N+M-1);
  w(q) = u'*ek{q}(1:N,1:M)*v;
end

% But if the basis functions aren't normalized, then we need to normalize the coefficients
w = w./self.norm(self.range(length(w))).^2;
