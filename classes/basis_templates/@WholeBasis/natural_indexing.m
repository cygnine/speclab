function subs = natural_indexing(self, subs)
% natural_indexing -- Whole number indexing
%
% inds = natural_indexing(self, subs)
%
%     The input array subs as entries in the natural numbers 1, 2, 3, ... and
%     the output "inds" replaces each input with its "whole number" index:
%
%      0 ------> 1
%      1 ------> 2
%      2 ------> 3
%      3 ------> 4
%      4 ------> 5
%           .
%           .
%           .

subs = subs + 1;
