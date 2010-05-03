function[x,w] = map_quadrature_rule(x,w,I1,I2)
% map_quadrature_rule -- Maps a 1D quadrature rule from one interval to another
%
% [x,w] = map_quadrature_rule(x,w,I1,I2)
%
%     Maps a quadrature rule (x,w) defined on 2-vector I1 to a new quadrature
%     rule on the 2-vector I2 via an affine transformation. The Jacobian is
%     applied to the new weights; no checking is performed to ensure that the
%     input nodes x are inside I1. x and w must be column vectors.
%
%     This function is `vectorized' (i.e. support multidimensional rules) in the
%     following way: if x and w are N x d matrices, and I1 and I2 are 2 x d
%     matrices, then this function performs the same 1D operation for each
%     column.

persistent spdiag
if isempty(spdiag)
  from labtools import spdiag
end

d = size(x,2);
if d==1
  scale = diff(I2)/diff(I1);

  x = (x-I1(1))*scale + I2(1);
  w = w*scale;
else
  scale = diff(I2, 1, 1).';
  x = (x-ones(size(x))*spdiag(I1(1,:).'))*spdiag(scale) + ...
      ones(size(x))*spdiag(I2(1,:).');

  if size(w,2)==1
    w = w.*prod(scale);
  else
    w = w*spdiag(scale);
  end
end
