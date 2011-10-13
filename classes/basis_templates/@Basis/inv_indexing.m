function[n] = inv_indexing(self, n_array)
% indexing -- general indexing method for class Basis
%
% [n] = inv_indexing(self, n_array)
%
%     Given a set of 'internal' (backend) indices n_array, this function
%     translates these indices into 'internal' (backend) indices for use in
%     computations. The basic idea is that n_array is (a) fed through
%     self.internal_indexing.inv and then (b) fed through self.user_indexing.
%
%     n_array should be an array where each row is an 'index'.

%self.internal_indexing.validate(n_array);
%self.internal_indexing.image.validate(n_array);

n = self.internal_indexing.inv(n_array);
% Just double-check
assert(all(n) > 0, 'Given indices do not match instance indexing rule');
n = self.user_indexing(n);

%self.user_indexing.image.validate(n);
%
%n_array = self.user_indexing.inv(n);
%numeln = numel(n_array);
%nsize = size(n_array);
%
%n_array = n_array(:);
%
%% Just double-check
%assert(all(n_array) > 0, 'Given indices do not match instance indexing rule');
%
%n_array = self.internal_indexing(n_array(:));
