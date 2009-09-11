function[y] = newton_evaluate(x,c,z)
% [Y] = NEWTON_EVALUATE(X,C,Z)
%
%     Uses Horner's method of evaluation to compute the interpolant at
%     the points Z given the modal coefficients C corresponding to a Newton
%     polynomial basis at nodes X.
%     length(C) = n
%     length(X) = n-1 (any additional nodal locations are ignored)

n = length(c);
y = c(n)*ones(size(z));
for q = (n-1):(-1):1
  y = y.*(z-x(q)) + c(q);
end
