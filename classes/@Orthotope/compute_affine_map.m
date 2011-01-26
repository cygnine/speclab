function[m] = compute_affine_map(self, varargin)
% compute_affine_map -- Computes affine map between two orthotopes
%
% m = compute_affine_map(self, other)
% m = compute_affine_map(self, {'key', value, 'Orthotope', specification})
%
%     Computes the affine map that connects two orthotopes,
%     dimension-by-dimension. Both orthotopes must have the same dimension, and
%     for each dimension, the lengths of the domain and target must both be
%     either finite or infinite. This generates the map m: self ---> other.

if not(isa(varargin{1}, 'Orthotope'))
  other = Orthotope(varargin{:});
else
  other = varargin{1};
end

if not(self.dimension == other.dimension)
  error('Both orthotopes must be of the same dimension');
end

A = zeros(self.dimension);
b = zeros([self.dimension 1]);

for q = 1:self.dimension
  flag0 = isfinite(self.slices{q}.length);
  flag1 = isfinite(other.slices{q}.length);
  if xor(flag0, flag1)
    error('There is no finite affine map connecting infinite and finite intervals')
  elseif not(flag1)
    % Connecting infinite intervals...don't really know what to do, just set
    % scale=1, define shift appropriately.
    A(q,q) = 1;
    b(q) = other.slices{q}.centroid - self.slices{q}.centroid;
  else
    % Connecting finite intervals, we're in business
    A(q,q) = other.slices{q}.length/self.slices{q}.length;
    b(q) = other.slices{q}.centroid - A(q,q)*self.slices{q}.centroid;
  end
end

m = AffineMap(A, b);
