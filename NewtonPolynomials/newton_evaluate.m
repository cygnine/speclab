function[y] = newton_evaluate(z,c,x)
% [Y] = NEWTON_EVALUATE(Z,C,X)
%
%     Uses Horner's method of evaluation to compute the interpolant at
%     the points z given the modal coefficients c corresponding to a Newton
%     polynomial basis at nodes x.
%     length(c) = n
%     length(x) = n-1 (any additional nodal locations are ignored)

n = length(c);
y = c(n)*ones(size(z));
for q = (n-1):(-1):1
  y = y.*(z-x(q)) + c(q);
end
