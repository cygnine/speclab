classdef WholeMultiIndices < IndexSet
% WholeMultiIndices -- The set of whole-number-valued multi-indices of a given dimension
%
% self = WholeMultiIndices(dim)
%
%     The set of multi-indices in dim-dimensional Euclidean space with
%     values in WholeNumbers for each dimension. This is a class with a
%     singleton instance for each dim.
%
%     An element of this set is a column vector with dim entries. An array
%     is a collection of these elements if it has dim rows.
%
% WholeMultiIndices Properties:
%   descriptor - A human-readable description: 'Whole-number-valued dim-dimenaional multi-indices'
% WholeMultiIndices Methods:
%   validate - Determines is an input is from the index set

  properties(SetAccess=private)
    descriptor
    dim
  end
  methods(Static)
    function self = instance(dim)
      persistent objs
      if nargin < 1
        dim = 1;
      end
      if isempty(objs)
        objs = {};
      end
      if (length(objs) < dim) || isempty(objs{dim})
        objs{dim} = WholeMultiIndices(dim);
      end
      self = objs{dim};
    end
  end
  methods(Access=private)
    function self = WholeMultiIndices(dim)
      self.dim = dim;
      self.descriptor = ['Whole-number-valued, ' num2str(dim) '-dimensional multi-indices'];
    end
  end
  methods
    validate(self,inp)
  end
end
