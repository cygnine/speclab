function[m] = compute_affine_map(self, varargin)
% compute_affine_map -- Computes affine map between two intervals
%
% m = compute_affine_map(self, other)
% m = compute_affine_map(self, {'key', value, 'Interval1D', specification})
%
%     Computes the affine map that connects two one-dimensional intervals. The
%     lengths of the domain and target must both be either finite or infinite.
%     This generates the map m: self ---> other.

if not(isa(varargin{1}, 'Interval1D'))
  other = Interval1D(varargin{:});
else
  other = varargin{1};
end

%m = inv(self.map_to_standard_interval);
%m = m.compose(other.map_to_standard_interval);

m = inv(other.map_to_standard_interval);
m = m.compose(self.map_to_standard_interval);
