function inds = index_natural_to_linear(self, inds)
% index_natural_to_linear -- 
%
% inds = index_natural_to_linear(self, inds)
%
%     The input inds is an array of entries with whole-basis indicies: 0, 1, 2,
%     etc. The output is the mapping from this set to the set of natural
%     numbers. In particular:
%
%      1 ------> 0
%      2 ------> 1
%      3 ------> 2
%      4 ------> 3
%      5 ------> 4
%           .
%           .
%           .

inds = inds - 1;
