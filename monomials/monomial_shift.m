function[c] = monomial_shift(c, shift)
% monomial_shift -- Returns coefficients of shifted polynomial
%
% c = monomial_shift(c,shift)
%
%     The input monomial is defined by N coefficients c: f(x) = c(1) + c(2)*x + ...
%     c(N)*x^(N-1). The returned coefficients are the monomial expansion of the
%     new function f(x-shift).

N = length(c);
monom_exp = zeros([N 1]);  
monom_exp(1) = 1;

% monom_exp represents the expansion of (x-shift)^q in monomials. Right now, q=0

oldc = c;
c = 0*c;

for q = 0:(N-1);
  c = c + oldc(q+1)*monom_exp;

  % Get ready for next q:
  monom_exp = -shift*monom_exp + [0; monom_exp(1:(N-1))];
end
