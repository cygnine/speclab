% MATLAB File : newton_evaluate.m
% [y] = newton_evaluate(z,c,x)
%
% * Creation Date : 2009-06-03
%
% * Last Modified : Fri 12 Jun 2009 03:21:07 PM EDT
%
% * Created By : Akil Narayan
%
% * Purpose : Uses Horner's method of evaluation to compute the interpolant at
%   the points z given the modal coefficients c corresponding to a Newton
%   polynomial basis at nodes x.
%   length(c) = n
%   length(x) = n-1 (any additional nodal locations are ignored)

function[y] = newton_evaluate(z,c,x)

n = length(c);
y = c(n)*ones(size(z));
for q = (n-1):(-1):1
  y = y.*(z-x(q)) + c(q);
end
