function[inds] = index_linear_to_natural(self, inds)
% index_linear_to_natural -- Converts linear indexing to the basis indexing
%
% inds = index_linear_to_natural(self, inds)
%
%     Given an array of linear indices (1, 2, 3, ...), this function coverts
%     these indices to the corresponding basis indexing (in this case, integer
%     indexing). Concretely:
%
%      1 ------>  0
%      2 ------> -1 
%      3 ------>  1
%      4 ------> -2 
%      5 ------>  2
%      6 ------> -3 
%      7 ------>  3
%      8 ------> -4 
%      9 ------>  4
%          .
%          .
%          .

inds = floor(inds/2).*(-1).^(inds+1);
