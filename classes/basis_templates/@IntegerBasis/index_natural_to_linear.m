function[inds] = index_natural_to_linear(self, inds)
% index_natural_to_linear -- Converts basis indexing to linear indexing
%
% inds = index_natural_to_linear(self, inds)
%
%     Given an array of basis-natural indices (0, -1, 1, -2, 2, ...) this
%     function coverts these indices to the corresponding lienar indexing.
%     Concretely:
%
%       0 ------> 1
%      -1 ------> 2
%       1 ------> 3
%      -2 ------> 4
%       2 ------> 5
%      -3 ------> 6
%       3 ------> 7
%      -4 ------> 8
%       4 ------> 9
%          .
%          .
%          .

inds = (sign(inds)>=0) + 2*abs(inds);
