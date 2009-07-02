% MATLAB File : newton_derivative_evaluate.m
% [f] = newton_derivative_evaluate(x,c)
%
% * Creation Date : 2009-06-03
%
% * Last Modified : Fri 12 Jun 2009 03:20:39 PM EDT
%
% * Created By : Akil Narayan
%
% * Purpose : For nodal locations x and modal coefficients c, calculates the
%   value of the derivative of the Newton interpolant at x(1). The Horner
%   decomposition of the interpolant makes this task very easy. 
%   If c has multiple columns (say C of them), we assume that there are multiple
%   interpolants whose derivatives must be evaluated simultaneously.
%   length(c) = n
%   length(x) = n-1  (any additional nodal locations are ignored)

function[f] = newton_derivative_evaluate(x,c)

[n,C] = size(c);
if and(n==1,C>1)  % I don't think you're calling this for derivatives of
  n = C;          % constants; I'm flipping some vectors around.
  C = 1;
  x = x';
  y = y';
end

f = c(n,:);
z = x(1,:);
for q = (n-1):(-1):2
  f = f.*(z-x(q,:)) + c(q,:);
end
