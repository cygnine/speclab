function[m] = compute_affine_map(self, varargin)
% compute_affine_map -- Computes affine map between two balls
%
% m = compute_affine_map(self, other)
% m = compute_affine_map(self, {'key', value, 'EuclideanBall', specification})
%
%     Computes the affine map that connects two Euclidean balls of the same
%     dimension. The radii of the domain and target must both be either finite
%     or infinite.  This generates the AffineMap instance map m: self --->
%     other.

if not(isa(varargin{1}, 'EuclideanBall'))
  other = EuclideanBall(varargin{:});
else
  other = varargin{1};
end

m = inv(other.map_to_standard_interval);
m = m.compose(self.map_to_standard_interval);
