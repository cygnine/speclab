function[self] = fullroots_from_modified_roots(self)
% fullroots_from_modified_roots -- Computes full set of roots
%
% roots = fullroots_from_modified_roots(self)
%
%     Saves a bit of typing by modularizing code that computes the full set of
%     roots from the exterior, conjugate, and quadratic categorization.

roots = zeros([self.num_exterior + self.num_conjugate + self.num_quadratic 1]);

rcount = 1;
for q = 1:self.num_exterior
  roots(rcount) = self.exterior_roots(q);
  rcount = rcount + 1;
end
for q = 1:self.num_conjugate
  roots(rcount) = self.conjugate_roots(q);
  roots(rcount+1) = conj(self.conjugate_roots(q));
  rcount = rcount + 2;
end
for q = 1:self.num_quadratic
  roots([rcount rcount+1]) = self.quadratic_roots(q);
  rcount = rcount + 2;
end

self.roots = roots;
