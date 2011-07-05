function[n_array, nsize, numeln] = translate_user_indexing(self,n)
% translate_user_indexing -- Translates user indexing into internal indexing
%
% [n_array, nsize, numeln] = translate_user_indexing(self,n)
%
%     Given an array n of user-type indices, this function translates them into
%     basis-internal indices, given in the output n_array. Since the internal
%     indexing for OrthogonalPolynomialBasis instances is ZeroBasedIndexing, we
%     can also return the size of the set of indices (nsize) and the number of
%     total index'ings (numeln). The output n_array is a column vector
%     containing the indices. 
%     (To reshape to original size do e.g. reshape(n_array, nsoze))

n_array = self.internal_indexing(self.user_indexing.inv(n));
nsize = size(n_array);
n_array = n_array(:);
numeln = length(n_array);
