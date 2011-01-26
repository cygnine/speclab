function[m] = inv(self)
% inv -- Inverse method for AffineMap
%
% m = inv(self)
%
%     Computes the inverse of the affine map, if the inverse exists.

if self.m ~= self.n
  error('Inverse does not exist');
elseif cond(self.A) > 1e20
  warning('Map is ill-conditioned, inverse may not be accurate');
end

Ainv = inv(self.A);
m = AffineMap(Ainv, -Ainv*self.b);
