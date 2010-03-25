function[x] = affine_transformation(x, scale, shift)
% affine_transformation -- Performs an affine transformation
%
% y = affine_transformation(x, scale, shift)
%
%     Performs an affine transformation in d dimensions, where d = size(x,2). x
%     is a collection of n points in d dimensions, it is of size (n, d). For q =
%     1, 2,..., d, x(:,q) is transformed under the affine transformation defined
%     by shift(q) and scale(q), which are both vectors of length d.

persistent spdiag
if isempty(spdiag)
  from labtools import spdiag
end

d = size(x,2);
if (length(scale)==1) & (length(shift)==1)
  x = scale*x+shift;
elseif length(scale)==1
  for q = 1:d
    x(:,q) = x(:,q) + shift(q);
  end
  x = x*scale;
elseif length(shift)==1
  x = (x + shift)*spdiag(scale);
else
  for q = 1:d
    x(:,q) = scale(q)*(x(:,q) + shift(q));
  end
end
