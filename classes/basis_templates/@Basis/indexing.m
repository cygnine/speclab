function[n_array, nsize, numeln] = indexing(self, n)
% indexing -- general indexing method for class Basis
%
% [n_array, nsize, numeln] = indexing(self, n)
%
%     Given a set of indices n, this function translates these indices into
%     'internal' (backend) indices for use in computations. The basic idea is
%     that n is (a) fed through self.user_indexing.inv and then (b) fed through
%     self.internal_indexing.
%
%     n_array contains the resulting (generally array) of indices. nsize is the
%     size of the array of input indices n, and numeln is the total number of
%     indices that n represents. The computation of nsize and n_array depends on
%     the defined indexing of self.

self.user_indexing.image.validate(n);

n_array = self.user_indexing.inv(n);
numeln = numel(n_array);
nsize = size(n_array);

n_array = n_array(:);

% Just double-check
assert(all(n_array) > 0, 'Given indices do not match instance indexing rule');

n_array = self.internal_indexing(n_array(:));
