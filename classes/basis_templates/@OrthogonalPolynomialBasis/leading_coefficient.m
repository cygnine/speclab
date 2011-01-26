function[k] = leading_coefficient(self, n)
% leading_coefficient -- Leading coefficient for polynomials
%
% k = leading_coefficient(self, n)
%
%     Returns the leading coefficients k of the degree-n orthogonal polynomials.

% Like self.l2_norm, this is purely a function of how things are scaled.
nsize = size(n);
n = n(:);
[a,b] = self.recurrence(0:(max(n)));
b = cumprod(sqrt(b));
k = 1./b(n+1);

% The above are the leading coefficients for the unmapped normalized polynomials
k = k.*self.scale_functions(ones(size(n)), n);
k = k.*(self.map_to_domain.A).^n;

k = reshape(k, nsize);
